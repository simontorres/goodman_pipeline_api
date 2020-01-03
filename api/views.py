from asgiref.sync import async_to_sync
from astropy.io import fits
from ccdproc import CCDData
import copy
import ccdproc
import matplotlib.pyplot as plt
import glob
import json
import logging
import re
import os
import numpy as np
import requests
import time
from channels.layers import get_channel_layer
from django.db import IntegrityError
from django.core import serializers

from ccdproc import ImageFileCollection

import goodman_pipeline

from goodman_pipeline.spectroscopy import WavelengthCalibration

from goodman_pipeline.core import (astroscrappy_lacosmic,
                                   call_cosmic_rejection,
                                   create_master_bias,
                                   create_master_flats,
                                   define_trim_section,
                                   extraction,
                                   get_slit_trim_section,
                                   get_overscan_region,
                                   identify_targets,
                                   image_overscan,
                                   image_trim,
                                   normalize_master_flat,
                                   read_fits,
                                   setup_logging,
                                   trace_targets,
                                   write_fits
                                   )

from goodman_pipeline.core import NoMatchFound

# from .models import ApiSettings
# from .models import ObservedFile
# from .models import ApiSessions, Calibrations

from rest_framework import status
from rest_framework.decorators import api_view
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework.permissions import IsAuthenticated
from .version import __version__

KEYWORDS = ['OBSTYPE',
            'INSTCONF',
            'CAM_TARG',
            'GRT_TARG',
            'FILTER',
            'FILTER2',
            'GRATING',
            'SLIT',
            'WAVMODE',
            'RDNOISE',
            'GAIN',
            'ROI']

setup_logging(debug=True)

log = logging.getLogger(__name__)


def add_trace_info(ccd, trace_info):
    last_keyword = None
    for keyword in trace_info:
        if last_keyword is None:
            ccd.header.set(keyword,
                           value=trace_info[keyword][0],
                           comment=trace_info[keyword][1])
        else:
            ccd.header.set(keyword,
                           value=trace_info[keyword][0],
                           comment=trace_info[keyword][1],
                           after=last_keyword)
        last_keyword = keyword
    return ccd


def auto_discover_files(full_path, file_type):

    if os.path.exists(full_path):
        if not os.path.isabs(full_path):
            full_path = os.path.abspath(full_path)

        ifc = ImageFileCollection(location=full_path, keywords=KEYWORDS)
        pd_ifc = ifc.summary.to_pandas()

        if file_type == 'BIAS':
            df_bias = pd_ifc[pd_ifc['OBSTYPE'] == 'BIAS']
            groups = df_bias.groupby(['OBSTYPE',
                                      'GAIN',
                                      'RDNOISE',
                                      'ROI']).size().reset_index().rename(
                columns={0: 'count'})
            bias_list = []
            for i in groups.index:
                bias_group = df_bias[
                    ((df_bias['OBSTYPE'] == groups.iloc[i]['OBSTYPE']) &
                     (df_bias['GAIN'] == groups.iloc[i]['GAIN']) &
                     (df_bias['RDNOISE'] == groups.iloc[i]['RDNOISE']) &
                     (df_bias['ROI'] == groups.iloc[i]['ROI']))]
                bias_list.append(bias_group.file.tolist())
            return bias_list

        if file_type == 'FLAT':
            df_flat = pd_ifc[pd_ifc['OBSTYPE'] == 'LAMPFLAT']
            groups = df_flat.groupby(['OBSTYPE',
                                      'FILTER',
                                      'FILTER2',
                                      'GRATING',
                                      'SLIT',
                                      'GAIN',
                                      'CAM_TARG',
                                      'GRT_TARG',
                                      'RDNOISE',
                                      'ROI']).size().reset_index().rename(
                columns={0: 'count'})
            flat_list = []
            for i in groups.index:
                flat_group = df_flat[
                    ((df_flat['OBSTYPE'] == groups.iloc[i]['OBSTYPE']) &
                     (df_flat['FILTER'] == groups.iloc[i]['FILTER']) &
                     (df_flat['FILTER2'] == groups.iloc[i]['FILTER2']) &
                     (df_flat['GRATING'] == groups.iloc[i]['GRATING']) &
                     (df_flat['SLIT'] == groups.iloc[i]['SLIT']) &
                     (df_flat['CAM_TARG'] == groups.iloc[i]['CAM_TARG']) &
                     (df_flat['GRT_TARG'] == groups.iloc[i]['GRT_TARG']) &
                     (df_flat['GAIN'] == groups.iloc[i]['GAIN']) &
                     (df_flat['RDNOISE'] == groups.iloc[i]['RDNOISE']) &
                     (df_flat['ROI'] == groups.iloc[i]['ROI']))]
                flat_list.append(flat_group.file.tolist())
            return flat_list

    else:
        raise FileNotFoundError('File {} does not exist'.format(full_path))


# def find_matching_calibration(full_path, calibration_type, search_path):
#     log.debug("Full Path: {} Calibration Type: {} Search Path: {}".format(
#         full_path,
#         calibration_type,
#         search_path))
#     if os.path.isdir(os.path.dirname(full_path)) and os.path.exists(full_path):
#         if not os.path.isabs(full_path):
#             full_path = os.path.abspath(full_path)
#
#         header = fits.getheader(full_path)
#         ifc = ImageFileCollection(location=os.path.abspath(search_path),
#                                   keywords=KEYWORDS)
#
#         pd_ifc = ifc.summary.to_pandas()
#
#         if calibration_type == 'BIAS':
#             filtered = pd_ifc[
#                 ((pd_ifc['OBSTYPE'] == calibration_type) &
#                  (pd_ifc['GAIN'] == header['GAIN']) &
#                  (pd_ifc['RDNOISE'] == header['RDNOISE']) &
#                  (pd_ifc['ROI'] == header['ROI'])
#                  )]
#             bias_list = filtered.file.tolist()
#             print(bias_list)
#
#
#             if len(bias_list) == 1:
#                 return bias_list[0]
#             elif len(bias_list) > 1:
#                 return bias_list[-1]
#
#             else:
#                 raise FileNotFoundError("Master bias not found for {}".format(
#                     os.path.basename(full_path)))
#
#         elif calibration_type == 'FLAT':
#             filtered = pd_ifc[
#                 (((pd_ifc['OBSTYPE'] == 'LAMPFLAT') |
#                   (pd_ifc['OBSTYPE'] == 'FLAT')) &
#                  (pd_ifc['FILTER'] == header['FILTER']) &
#                  (pd_ifc['FILTER2'] == header['FILTER2']) &
#                  (pd_ifc['GRATING'] == header['GRATING']) &
#                  (pd_ifc['SLIT'] == header['SLIT']) &
#                  (pd_ifc['CAM_TARG'] == header['CAM_TARG']) &
#                  (pd_ifc['GRT_TARG'] == header['GRT_TARG']) &
#                  (pd_ifc['GAIN'] == header['GAIN']) &
#                  (pd_ifc['RDNOISE'] == header['RDNOISE']) &
#                  (pd_ifc['ROI'] == header['ROI'])
#                  )]
#             flat_list = filtered.file.tolist()
#
#             print("FLAT LIOST")
#             print(flat_list)
#             if len(flat_list) == 1:
#                 return flat_list[0]
#             elif len(flat_list) == 2:
#                 for master in flat_list:
#                     if 'norm_master' in master:
#                         return master
#             else:
#                 raise FileNotFoundError("Master flat not found for {}".format(
#                     os.path.basename(full_path)))
#
#     else:
#         raise FileNotFoundError("Reference File {} does not exist".formata(
#             os.path.basename(full_path)))


class ApiView(APIView):
    def get(self, request):
        content = {'detail': 'Data Reduction Goodman Pipeline API'}
        return Response(content)


class ApiCalibrationsBiasView(APIView):

    response_payload = {
        'api_version': __version__,
        'pipeline_version': goodman_pipeline.__version__,
        'request_data': '',
        'error': '',
        'master_bias': ''}
    #
    # def get(self, request):
    #     response_payload = copy.deepcopy(self.response_payload)
    #     file_name = request.GET.get('file_name', '')
    #     raw_data_path = request.GET.get('raw_data_path')
    #     reduced_data_path = request.GET.get('reduced_data_path')
    #
    #     response_payload['request_data'] = {
    #         'file_name': file_name,
    #         'raw_data_path': raw_data_path,
    #         'reduced_data_path': reduced_data_path}
    #
    #     if file_name != '':
    #         try:
    #             raw_full_path = os.path.join(raw_data_path, file_name)
    #             reduced_full_path = os.path.join(reduced_data_path,
    #                                              file_name)
    #             file_full_path = next(
    #                 path for path in [raw_full_path, reduced_full_path] if
    #                 os.path.isfile(path))
    #
    #             master_bias_file = find_matching_calibration(
    #                 full_path=file_full_path,
    #                 calibration_type='BIAS',
    #                 search_path=reduced_data_path)
    #             print(master_bias_file)
    #             response_payload['master_bias'] = master_bias_file
    #
    #             return Response(response_payload)
    #         except FileNotFoundError as error:
    #             response_payload['error'] = str(error)
    #             return Response(response_payload,
    #                             status=status.HTTP_404_NOT_FOUND)
    #         except StopIteration:
    #             response_payload['error'] = "File {} not found".format(
    #                 file_name)
    #             return Response(response_payload,
    #                             status=status.HTTP_404_NOT_FOUND)
    #     else:
    #
    #         all_master_bias = auto_discover_files(full_path=reduced_data_path,
    #                                               file_type='BIAS')
    #         if len(all_master_bias) == 1:
    #             response_payload['master_bias'] = all_master_bias
    #             return Response(response_payload)
    #         else:
    #             response_payload['error'] = "Unable to find master bias files"
    #             return Response(response_payload,
    #                             status=status.HTTP_404_NOT_FOUND)

    def post(self, request):
        response_payload = copy.deepcopy(self.response_payload)
        auto_discover = request.data.get('auto_discover', 'true')
        raw_data_path = request.data.get('raw_data_path')
        reduced_data_path = request.data.get('reduced_data_path')
        file_list = request.data.get('file_list', [])
        technique = request.data.get('technique')

        response_payload['request_data'] = {
            'auto_discover': auto_discover,
            'raw_data_path': raw_data_path,
            'reduced_data_path': reduced_data_path,
            'file_list': file_list,
            'technique': technique}

        if auto_discover == 'true':
            try:
                file_list = auto_discover_files(full_path=raw_data_path,
                                                file_type='BIAS')
            except FileNotFoundError as error:
                response_payload['error'] = str(error)
                return Response(response_payload,
                                status=status.HTTP_404_NOT_FOUND)
        if file_list != []:
            all_master_bias = []
            for sub_list in file_list:
                if isinstance(sub_list, list):
                    print("is a list ")
                    start_time = time.time()
                    master_bias, master_bias_name = create_master_bias(
                        bias_files=sub_list,
                        raw_data=raw_data_path,
                        reduced_data=reduced_data_path,
                        technique=technique)
                    elapsed = time.time() - start_time

                    all_master_bias.append({
                        'parent_file': os.path.basename(master_bias_name),
                        'file_name': os.path.basename(master_bias_name),
                        'full_path': os.path.join(reduced_data_path,
                                                  master_bias_name),
                        'file_list': sub_list,
                        'elapsed_time': elapsed})
                else:
                    print("not a list")
                    print(dir(sub_list))

            response_payload['master_bias'] = all_master_bias

            return Response(response_payload)
        else:
            response_payload['error'] = 'List of files is empty'
            return Response(response_payload)


class ApiCalibrationsFlatsView(APIView):
    response_payload = {
        'api_version': __version__,
        'pipeline_version': goodman_pipeline.__version__,
        'request_data': '',
        'error': '',
        'master_flats': ''}

    # def get(self, request):
    #     response_payload = copy.deepcopy(self.response_payload)
    #
    #     file_full_path = request.GET.get('file_full_path', '')
    #     reduced_data_path = request.GET.get('reduced_data_path', '')
    #     master_bias_file = request.GET.get('master_bias_file')
    #     response_payload['request_data'] = {
    #         'file_full_path': file_full_path,
    #         'reduced_data_path': reduced_data_path}
    #
    #     if file_full_path != '':
    #         try:
    #             raw_full_path = os.path.join(raw_data_path, file_full_path)
    #             reduced_full_path = os.path.join(reduced_data_path, file_full_path)
    #             full_path = next(path for path in [raw_full_path, reduced_full_path] if os.path.isfile(path))
    #
    #             master_flat_file = find_matching_calibration(
    #                 full_path=full_path,
    #                 calibration_type='FLAT',
    #                 search_path=reduced_data_path)
    #
    #             response_payload['master_flats'] = master_flat_file
    #
    #             return Response(response_payload)
    #         except FileNotFoundError as error:
    #             response_payload['error'] = str(error)
    #             return Response(response_payload,
    #                             status=status.HTTP_404_NOT_FOUND)
    #         except StopIteration:
    #             response_payload['error'] = 'File {} does not exist' \
    #                                              ''.format(file_full_path)
    #             return Response(response_payload)
    #
    #     else:
    #         all_flats = auto_discover_files(full_path=reduced_data_path,
    #                                         file_type='FLAT')
    #         response_payload['master_flats'] = all_flats
    #         return Response(response_payload)

    def post(self, request):
        response_payload = copy.deepcopy(self.response_payload)
        auto_discover = request.data.get('auto_discover')
        raw_data_path = request.data.get('raw_data_path')
        reduced_data_path = request.data.get('reduced_data_path')
        overscan_region = request.data.get('overscan_region')
        trim_section = request.data.get('trim_section')
        master_bias_file = request.data.get('master_bias_file')
        file_list = request.data.get('file_list')
        technique = request.data.get('technique')
        settings = request.data.get('settings')

        flat_normalization = settings['flat_normalization']
        flat_normalization_order = settings['flat_normalization_order']
        saturation_threshold = settings['saturation_threshold']
        ignore_bias = settings['ignore_bias']

        response_payload['request_data'] = {
            'auto_discover': auto_discover,
            'raw_data_path': raw_data_path,
            'reduced_data_path': reduced_data_path,
            'overscan_region' : overscan_region,
            'trim_section': trim_section,
            'master_bias_file': master_bias_file,
            'file_list': file_list,
            'technique': technique,
            'settings': settings
        }

        log.debug("Master Bias File: {}".format(master_bias_file))
        if not ignore_bias:
            if master_bias_file is None:
                response_payload['error'] = 'Must provide a master_bias_' \
                                            'file or set ignore_bias to False'
                log.error(response_payload['error'])
                return Response(response_payload,
                                status=status.HTTP_400_BAD_REQUEST)
            elif not os.path.isfile(master_bias_file) or not os.path.exists(master_bias_file):
                response_payload['error'] = "Unable to locate master bias file "\
                                            "{}".format(master_bias_file)
                log.error(response_payload['error'])
                return Response(response_payload,
                                status=status.HTTP_404_NOT_FOUND)
        else:
            log.warning('Ignore bias set to {}'.format(ignore_bias))
            master_bias_file = ''
        if auto_discover:
            file_list = auto_discover_files(full_path=raw_data_path,
                                            file_type='FLAT')
        all_flats = []
        for sub_list in file_list:
            sample_file = os.path.join(raw_data_path,
                                       sub_list[0])

            if overscan_region == "":
                overscan_region = get_overscan_region(sample_image=sample_file,
                                                      technique=technique)
                response_payload['overscan_region'] = overscan_region
            if trim_section == "":
                trim_section = define_trim_section(sample_image=sample_file,
                                                   technique=technique)
                response_payload['trim_section'] = trim_section

            file_header = fits.getheader(sample_file)

            if technique == 'Spectroscopy':
                new_master_flat_name = "master_flat_{}_{}.fits".format(
                    re.sub(' ', '_', file_header['WAVMODE']),
                    re.sub('[<> ]', '', file_header['FILTER2']))
            else:
                new_master_flat_name = "master_flat_{}_{}.fits".format(
                    re.sub(' ', '_', file_header['WAVMODE']),
                    re.sub('[<> ]', '', file_header['FILTER']))


            try:
                master_flat, master_flat_name = create_master_flats(
                    flat_files=sub_list,
                    raw_data=raw_data_path,
                    reduced_data=reduced_data_path,
                    technique=technique,
                    overscan_region=overscan_region,
                    trim_section=trim_section,
                    master_bias_name=master_bias_file,
                    new_master_flat_name=new_master_flat_name,
                    saturation=saturation_threshold,
                    ignore_bias=ignore_bias)

                if technique == 'Spectroscopy':
                    log.info("Normalizing {} flat by method: {}".format(
                        technique, flat_normalization))
                    norm_master,  norm_master_name = normalize_master_flat(
                        master=master_flat,
                        name=master_flat_name,
                        method=flat_normalization,
                        order=int(flat_normalization_order))
                else:
                    log.info("Normalizing {} flat by method: {}".format(
                        technique, 'mean'))
                    norm_master,  norm_master_name = normalize_master_flat(
                        master=master_flat,
                        name=master_flat_name,
                        method='mean',
                        order=int(flat_normalization_order))

                all_flats.append({
                    'parent_file': os.path.basename(master_flat_name),
                    'file_name': os.path.basename(master_flat_name),
                    'full_path': master_flat_name,
                    'normalized': False,
                    'file_list': sub_list
                    })
                all_flats.append({
                    'parent_file': os.path.basename(master_flat_name),
                    'file_name': os.path.basename(norm_master_name),
                    'full_path': norm_master_name,
                    'normalized': True,
                    'file_list': sub_list
                    })

                response_payload['master_flats'] = all_flats
                return Response(response_payload)
            except ValueError as error:
                response_payload['error'] = str(error)
                log.error(response_payload['error'])
            except TypeError as error:
                response_payload['error'] = str(error)
                log.error(response_payload['error'])

            finally:
                return Response(response_payload)





class ApiImageOverscanView(APIView):

    def get(self, request):
        _file = request.GET.get('file')
        raw_data_path = request.GET.get('raw_data_path')
        full_path = os.path.join(raw_data_path, _file)
        if os.path.isfile(full_path) and os.path.exists(full_path):
            overscan_region = get_overscan_region(sample_image=full_path,
                                                  technique='Spectroscopy')
            return Response({'api_version': __version__,
                             'pipeline_version': goodman_pipeline.__version__,'overscan_region': overscan_region})
        else:
            return Response({'api_version': __version__,
                             'pipeline_version': goodman_pipeline.__version__,'error': 'File {} does not exist'.format(_file)},
                            status=status.HTTP_404_NOT_FOUND)

    def post(self, request):
        _file = request.data.get('file')
        raw_data_path = request.data.get('raw_data_path')
        reduced_data_path = request.data.get('reduced_data_path')
        full_path = os.path.join(raw_data_path, _file)
        if os.path.isfile(full_path) and os.path.exists(full_path):
            ccd = read_fits(full_path=full_path, technique='Spectroscopy')

            overscan_region = get_overscan_region(sample_image=full_path,
                                                  technique='Spectroscopy')

            overscan_corrected = image_overscan(ccd=ccd,
                                                overscan_region=overscan_region,
                                                add_keyword=False)
            new_file_name = os.path.join(reduced_data_path, 'o_{}'.format(_file))

            write_fits(ccd=overscan_corrected,
                       full_path=new_file_name,
                       combined=False,
                       parent_file=_file,
                       overwrite=True)

            return Response({'api_version': __version__,
                             'pipeline_version': goodman_pipeline.__version__,'overscan_region': overscan_region,
                             'overscan_corrected': os.path.basename(new_file_name),
                             'original_file': _file})
        else:
            return Response({'api_version': __version__,
                             'pipeline_version': goodman_pipeline.__version__,'error': 'File {} does not exist'.format(_file)},
                            status=status.HTTP_404_NOT_FOUND)


class ApiImageTrimView(APIView):

    def get(self, request):
        _file = request.GET.get('file')
        path_type = request.GET.get('path_type', 'reduced')
        raw_data_path = request.GET.get('raw_data_path')
        reduced_data_path = request.GET.get('reduced_data_path')
        trim_type = request.GET.get('trim_type', 'trimsec')

        if path_type == 'reduced':
            full_path = os.path.join(reduced_data_path, _file)
        elif path_type == 'raw':
            full_path = os.path.join(raw_data_path, _file)
        else:
            return Response({'api_version': __version__,
                             'pipeline_version': goodman_pipeline.__version__,'error': "Unrecognized data location {}".format(path_type)})
        if os.path.isfile(full_path) and os.path.exists(full_path):
            if trim_type == 'trimsec':
                trim_section = define_trim_section(sample_image=full_path,
                                                   technique='Spectroscopy')
                return Response({'api_version': __version__,
                             'pipeline_version': goodman_pipeline.__version__,'trim_section': trim_section})
            elif trim_type == 'slit':
                if fits.getval(filename=full_path, keyword='OBSTYPE') == 'FLAT':
                    master_flat = CCDData.read(full_path, unit='adu')
                    slit_trim_section = get_slit_trim_section(
                        master_flat=master_flat)
                    return Response({'api_version': __version__,
                             'pipeline_version': goodman_pipeline.__version__,'slit_trim_section': slit_trim_section})
                else:
                    return Response({'api_version': __version__,
                             'pipeline_version': goodman_pipeline.__version__,
                        'error': "for obtaining slit trim section "
                                 "an image of OBSTYPE == FLAT is "
                                 "required. Parsed a {} type"
                                 "".format(fits.getval(full_path,
                                                       'OBSTYPE'))})
        else:
            return Response(
                {'error': 'File {} does not exist'.format(_file)},
                status=status.HTTP_404_NOT_FOUND)

    def post(self, request):
        _file = request.data.get('file')
        trim_section = request.data.get('trim_section')
        raw_data_path = request.data.get('raw_data_path')
        reduced_data_path = request.data.get('reduced_data_path')
        path_type = request.data.get('path_type')
        trim_type = request.data.get('trim_type', 'trimsec')



        if path_type == 'reduced':
            full_path = os.path.join(reduced_data_path, _file)
        elif path_type == 'raw':
            full_path = os.path.join(raw_data_path, _file)
        else:
            return Response({'api_version': __version__,
                             'pipeline_version': goodman_pipeline.__version__,'error': "Unrecognized data location {}".format(path_type)})

        if os.path.isfile(full_path) and os.path.exists(full_path):
            if trim_section == "":
                trim_section = define_trim_section(sample_image=full_path,
                                                   technique="Spectroscopy")

            ccd = read_fits(full_path=full_path, technique='Spectroscopy')
            trimmed_ccd = image_trim(ccd=ccd,
                                     trim_section=trim_section,
                                     trim_type=trim_type,
                                     add_keyword=False)

            new_file_name = os.path.join(reduced_data_path,
                                         't{}'.format(_file))

            write_fits(ccd=trimmed_ccd,
                       full_path=new_file_name,
                       combined=False,
                       parent_file=_file,
                       overwrite=True)

            return Response({'api_version': __version__,
                             'pipeline_version': goodman_pipeline.__version__,'trim_section': trim_section,
                             'trimmed_image': os.path.basename(
                                 new_file_name),
                             'trim_type': trim_type,
                             'parent_file': _file})
        else:
            return Response(
                {'error': 'File {} does not exist'.format(_file)},
                status=status.HTTP_404_NOT_FOUND)


class ApiImageBiasView(APIView):

    def get(self, request):
        return Response({'api_version': __version__,
                             'pipeline_version': goodman_pipeline.__version__,})

    def post(self, request):
        parent_file = request.data.get('file', '')
        if parent_file == '':
            return Response({'api_version': __version__,
                             'pipeline_version': goodman_pipeline.__version__,'error': 'Must provide a valid file name'})

        reduced_data_path = request.data.get('reduced_data_path')
        master_bias_file = request.data.get('master_bias', '')

        if os.path.exists(reduced_data_path) and os.path.isdir(reduced_data_path):
            full_path = os.path.join(reduced_data_path, parent_file)
            if not os.path.exists(full_path) or not os.path.isfile(full_path):
                return Response(
                    {'error': 'file {} does not exist'.format(parent_file)})
            if os.path.isabs(master_bias_file) and os.path.isfile(master_bias_file):
                master_bias = read_fits(full_path=master_bias_file,
                                        technique='Spectroscopy')

                ccd = read_fits(full_path=full_path,
                                technique='Spectroscopy')

                bias_corrected = ccdproc.subtract_bias(ccd=ccd,
                                                       master=master_bias,
                                                       add_keyword=False)
                bias_corrected_name = "z{}".format(parent_file)

                bias_corrected.header.set('GSP_PNAM', value=parent_file)
                bias_corrected.header.set('GSP_FNAM', value=bias_corrected_name)
                bias_corrected.header.set('GSP_BIAS', value=master_bias_file)

                bias_corrected.write(os.path.join(reduced_data_path,
                                                  bias_corrected_name),
                                     overwrite=True)

                return Response({'api_version': __version__,
                         'pipeline_version': goodman_pipeline.__version__,'master_bias': master_bias_file,
                                 'bias_corrected': bias_corrected_name})

            elif not os.path.exists(os.path.join(reduced_data_path,
                                                 master_bias_file)):
                return Response({'api_version': __version__,
                             'pipeline_version': goodman_pipeline.__version__,
                    'error': 'master bias file {} '
                             'does not exist'.format(master_bias_file)})
        else:
            return Response({'api_version': __version__,
                             'pipeline_version': goodman_pipeline.__version__,
                'error': 'folder {} does not exist'.format(
                    reduced_data_path)})


class ApiImageFlatView(APIView):

    def get(self, request):
        return Response({'api_version': __version__,
                             'pipeline_version': goodman_pipeline.__version__,'error': "Not implemented"})

    def post(self, request):
        parent_file = request.data.get('file', '')
        if parent_file == '':
            return Response({'api_version': __version__,
                             'pipeline_version': goodman_pipeline.__version__,'error': 'Must provide a valid file name'})

        reduced_data_path = request.data.get('reduced_data_path')
        master_flat_file = request.data.get('master_flat', '')

        if os.path.exists(reduced_data_path) and os.path.isdir(
                reduced_data_path):
            full_path = os.path.join(reduced_data_path, parent_file)
            if not os.path.exists(full_path) or not os.path.isfile(full_path):
                return Response(
                    {'error': 'file {} does not exist'.format(parent_file)})
            if master_flat_file == '':
                try:
                    master_flat_file = find_matching_calibration(
                        full_path=full_path, calibration_type='FLAT',
                        search_path=reduced_data_path)
                    master_flat = read_fits(
                        full_path=os.path.join(
                            reduced_data_path,
                            master_flat_file),
                        technique='Spectroscopy')

                    ccd = read_fits(full_path=full_path,
                                    technique='Spectroscopy')

                    flat_corrected = ccdproc.flat_correct(ccd=ccd,
                                                          flat=master_flat,
                                                          add_keyword=False)

                    flat_corrected_name = "f{}".format(parent_file)

                    flat_corrected.header.set('GSP_PNAM', value=parent_file)
                    flat_corrected.header.set('GSP_FNAM',
                                              value=flat_corrected_name)
                    flat_corrected.header.set('GSP_BIAS',
                                              value=master_flat_file)

                    flat_corrected.write(os.path.join(reduced_data_path,
                                                      flat_corrected_name),
                                         overwrite=True)

                    return Response({'api_version': __version__,
                             'pipeline_version': goodman_pipeline.__version__,'master_flat': master_flat_file,
                                     'flat_corrected': flat_corrected_name})
                except FileNotFoundError:
                    return Response({'api_version': __version__,
                             'pipeline_version': goodman_pipeline.__version__,
                        'error': "Unable to find a suitable master "
                                 "bias for file {}".format(parent_file)})
            elif not os.path.exists(os.path.join(reduced_data_path,
                                                 master_flat_file)):
                return Response({'api_version': __version__,
                                 'pipeline_version': goodman_pipeline.__version__,
                                 'error': 'master bias file {} '
                                          'does not exist'.format(master_flat_file)})
        else:
            return Response({'api_version': __version__,
                             'pipeline_version': goodman_pipeline.__version__,
                             'error': 'folder {} does not exist'.format(
                    reduced_data_path)})


class ApiImageCrejectView(APIView):

    def post(self, request):
        response_payload = {'api_version': __version__,
                            'pipeline_version': goodman_pipeline.__version__,
                            'error': '',
                            'request_data': {},
                            'crejected': ''}

        filename = request.data.get('filename', '')
        raw_data_path = request.data.get('raw_data_path')
        reduced_data_path = request.data.get('reduced_data_path')
        method = request.data.get('method', 'default')
        technique = request.data.get('technique', 'Spectroscopy')

        response_payload['request_data'] = {
            'filename': filename,
            'raw_data_path': raw_data_path,
            'reduced_data_path': reduced_data_path,
            'method': method,
            'technique': technique
        }

        if filename != '':
            raw_full_path = os.path.join(raw_data_path, filename)
            reduced_full_path = os.path.join(reduced_data_path, filename)

            full_path = next(path for path in [raw_full_path, reduced_full_path] if os.path.isfile(path))
            if not os.path.exists(full_path) or not os.path.isfile(full_path):
                response_payload['error'] = 'file {} does not exist'.format(filename)
                return Response(response_payload,
                                status=status.HTTP_404_NOT_FOUND)

            ccd = read_fits(full_path=full_path, technique=technique)

            ccd, _ = call_cosmic_rejection(
                ccd=ccd,
                image_name=filename,
                out_prefix='',
                red_path=reduced_data_path,
                keep_files=False,
                prefix='c',
                method=method,
                save=False)

            new_file_name = "c{}".format(filename)
            full_path = os.path.join(reduced_data_path, new_file_name)

            write_fits(ccd=ccd, full_path=full_path)
            response_payload['crejected'] = os.path.basename(full_path)

            return Response(response_payload)

        else:
            response_payload['error'] = 'Must provide a valid file name'
            return Response(response_payload,
                            status=status.HTTP_400_BAD_REQUEST)


class ApiImageReduceView(APIView):

    def post(self, request):
        response_payload = {'api_version': __version__,
                            'pipeline_version': goodman_pipeline.__version__,
                            'error': '',
                            'request_data': {}}

        file_full_path = request.data.get('file_full_path', '')
        reduced_data_path = request.data.get('reduced_data_path')
        master_bias_file = request.data.get('master_bias_file', '')
        master_flat_file = request.data.get('master_flat_file', '')
        technique = request.data.get('technique')
        settings = request.data.get('settings')

        cosmic_ray_rejection = settings['cosmic_ray_rejection']
        ignore_bias = settings['ignore_bias']
        ignore_flats = settings['ignore_flats']

        response_payload['request_data'] = {
            'file_full_path': file_full_path,
            'reduced_data_path': reduced_data_path,
            'master_bias_file': master_bias_file,
            'master_flat_file': master_flat_file,
            'technique': technique,
            'setttings': settings
        }
        if file_full_path == '' or not os.path.isabs(file_full_path):
            response_payload['error'] = 'please provide a full path to file'
            log.error(response_payload['error'])
            return Response(response_payload,
                            status=status.HTTP_400_BAD_REQUEST)

        if not os.path.exists(reduced_data_path) or not os.path.isdir(reduced_data_path):
            response_payload['error'] = 'reduced_data_path does not exist or ' \
                                        'is not a directory'
            log.error(response_payload['error'])
            return Response(response_payload,
                            status=status.HTTP_400_BAD_REQUEST)

        if not os.path.exists(file_full_path) or not os.path.isfile(file_full_path):
            response_payload['error'] = "{} does not exist in raw_data_path or " \
                                        "is not a file".format(file_full_path)
            log.error(response_payload['error'])
            return Response(response_payload,
                            status=status.HTTP_400_BAD_REQUEST)

        if not os.path.isabs(master_bias_file) or not os.path.isfile(master_bias_file):
            response_payload['error'] = "Must provide an absolute path for master_bias_file"
            log.error(response_payload['error'])
            return Response(response_payload,
                            status=status.HTTP_400_BAD_REQUEST)

        if not os.path.isabs(master_flat_file) or not os.path.isfile(master_flat_file):
            response_payload['error'] = "Must provide an absolute path for master_flat_file"
            log.error(response_payload['error'])
            return Response(response_payload,
                            status=status.HTTP_400_BAD_REQUEST)

        overscan_region = get_overscan_region(sample_image=file_full_path,
                                              technique=technique)

        trim_section = define_trim_section(sample_image=file_full_path,
                                           technique=technique)
        slit_trim_section = None
        if not ignore_flats:
            log.debug("Reading Master Flat File {}".format(master_flat_file))
            master_flat = read_fits(full_path=master_flat_file, technique=technique)

            master_flat = image_trim(ccd=master_flat,
                                     trim_section=trim_section,
                                     trim_type='trimsec',
                                     add_keyword=False)

            slit_trim_section = get_slit_trim_section(master_flat=master_flat)

            master_flat = image_trim(ccd=master_flat,
                                     trim_section=slit_trim_section,
                                     trim_type='slit',
                                     add_keyword=False)
        else:
            master_flat = None

        if not ignore_bias:
            log.debug("Reading Master Bias File {}".format(master_bias_file))
            master_bias = read_fits(full_path=master_bias_file,
                                    technique=technique)

            master_bias = image_trim(ccd=master_bias,
                                     trim_section=trim_section,
                                     trim_type='trimsec',
                                     add_keyword=False)

            if slit_trim_section is not None:
                master_bias = image_trim(ccd=master_bias,
                                         trim_section=slit_trim_section,
                                         trim_type='slit',
                                         add_keyword=False)
        else:
            master_bias = None

        prefix = '_'

        ccd = read_fits(full_path=file_full_path, technique=technique)

        log.info("Raw: Mean {}".format(np.mean(ccd.data)))

        if ignore_bias:
            ccd = image_overscan(ccd=ccd,
                                 overscan_region=overscan_region,
                                 add_keyword=False)
            prefix = 'o' + prefix
            log.info("Overscan: Mean {}".format(np.mean(ccd.data)))

        # save intermediate?

        ccd = image_trim(ccd=ccd,
                         trim_section=trim_section,
                         trim_type='trimsec')
        prefix = 't' + prefix

        log.info("Trim (trimsec): Mean {}".format(np.mean(ccd.data)))

        if slit_trim_section is not None:
            ccd = image_trim(ccd=ccd,
                             trim_section=slit_trim_section,
                             trim_type='slit')
            log.info("Trim (slit): Mean {}".format(np.mean(ccd.data)))

            prefix = 's' + prefix

        if not ignore_bias:
            print(ccd.data.shape, master_bias.data.shape)
            ccd = ccdproc.subtract_bias(ccd=ccd,
                                        master=master_bias,
                                        add_keyword=False)
            log.info("Bias: Mean {}".format(np.mean(ccd.data)))
            ccd.header.set('GSP_BIAS',
                           value=os.path.basename(master_bias_file))

            prefix = 'z' + prefix

        if not ignore_flats:
            ccd = ccdproc.flat_correct(ccd=ccd,
                                       flat=master_flat,
                                       add_keyword=False)
            log.info("Flat: Mean {}".format(np.mean(ccd.data)))
            ccd.header.set('GSP_FLAT', value=os.path.basename(master_flat_file))

            prefix = 'f' + prefix

        ccd, _ = call_cosmic_rejection(
            ccd=ccd,
            image_name=os.path.basename(file_full_path),
            out_prefix='staged_',
            red_path=reduced_data_path,
            keep_files=False,
            prefix='c',
            method=cosmic_ray_rejection,
            save=False)
        log.info("C Reject: Mean {}".format(np.mean(ccd.data)))
        prefix = 'c' + prefix

        new_name = prefix + os.path.basename(file_full_path)
        new_full_path = os.path.join(reduced_data_path, new_name)

        write_fits(ccd=ccd,
                   full_path=new_full_path,
                   parent_file=os.path.basename(file_full_path))

        # ccd.write(new_full_path, overwrite=True)

        response_payload['overscan_region'] = overscan_region
        response_payload['trim_section'] = trim_section
        response_payload['slit_trim_section'] = slit_trim_section
        response_payload['master_bias'] = master_bias_file
        response_payload['master_flat'] = master_flat_file
        response_payload['full_path'] = new_full_path
        response_payload['file_name'] = new_name

        return Response(response_payload)


class ApiSpectrumExtractView(APIView):

    def post(self, request):
        full_path = request.data.get('file_full_path', '')
        reduced_data_path = request.data.get('reduced_data_path',
                                             os.path.dirname(full_path))
        reference_lamps = request.data.get('reference_lamps', [])

        reduction_settings = request.data.get('settings')
        target_fit_model = reduction_settings['target_fit_model']
        background_threshold = reduction_settings['background_threshold']
        max_targets = reduction_settings['max_targets']
        extraction_type = reduction_settings['extraction_type']
        technique = 'Spectroscopy'

        response_payload = {'api_version': __version__,
                            'pipeline_version': goodman_pipeline.__version__,
                            'request_data':
                                {'file_full_path': full_path,
                                 'reduced_data_path': reduced_data_path,
                                 'reference_lamps': reference_lamps,
                                 'target_fit_model': target_fit_model,
                                 'background_threshold': background_threshold,
                                 'max_targets': max_targets,
                                 'extraction_type': extraction_type,
                                 'technique': technique
                                 },
                            'error': ''}

        if not os.path.isfile(full_path) or not os.path.exists(full_path):
            response_payload['error'] = 'Unable to find file {}'.format(full_path)
            log.critical(response_payload['error'])
            return Response(response_payload, status=status.HTTP_404_NOT_FOUND)

        if isinstance(reference_lamps, list):
            for _file in reference_lamps:
                if not os.path.isfile(_file) or not os.path.exists(_file):
                    response_payload['error'] = 'Unable to locate reference ' \
                                                'lamp file {}'.format(_file)
                    log.critical(response_payload['error'])
                    return Response(response_payload, status=status.HTTP_404_NOT_FOUND)
        else:
            response_payload['error'] = "reference_lamps must be a list"
            return Response(response_payload,
                            status=status.HTTP_400_BAD_REQUEST)

        if not os.path.exists(reduced_data_path) or not os.path.isdir(
                reduced_data_path):
            response_payload['error'] = 'data_path: {} does not exist or is ' \
                                        'not a directory'.format(reduced_data_path)
            log.critical(response_payload['error'])
            return Response(response_payload, status=status.HTTP_400_BAD_REQUEST)

        ccd = read_fits(full_path=full_path, technique=technique)
        try:
            targets = identify_targets(ccd=ccd,
                                       fit_model=target_fit_model,
                                       background_threshold=int(background_threshold),
                                       nfind=int(max_targets))

        except Exception as error:
            response_payload['error'] = '{}: {}'.format(
                error.__class__.__name__, str(error))
            log.error(response_payload['error'])
            return Response(response_payload, status=status.HTTP_500_INTERNAL_SERVER_ERROR)
        if len(targets) == 0:
            response_payload['error'] = 'Unable to identify targets in {}'.format(full_path)
            log.error(response_payload['error'])
            return Response(response_payload)

        traces = trace_targets(ccd=ccd,
                               target_list=targets,
                               sampling_step=5,
                               pol_deg=2,
                               nfwhm=5,
                               plots=False)

        i = 1
        all_targets = []
        for trace, profile, trace_info in traces:

            ccd = add_trace_info(ccd=ccd, trace_info=trace_info)

            target_name = "target_{}".format(i)
            extracted_ccd = extraction(ccd=ccd,
                                       target_trace=trace,
                                       spatial_profile=profile,
                                       extraction_name=extraction_type)
            new_file_name = "e{}".format(
                re.sub('.fits',
                       '_{}.fits'.format(int(trace.c0.value)),
                       os.path.basename(full_path)))

            new_full_path = os.path.join(reduced_data_path, new_file_name)
            write_fits(ccd=extracted_ccd,
                       full_path=new_full_path,
                       combined=False,
                       parent_file=full_path,
                       overwrite=True)

            this_target = {'target_center': trace.c0.value,
                           'file_full_path': new_full_path,
                           'reference_lamps':[],
                           'trace_rms_error': trace_info['GSP_TERR'][0]}



            for reference_lamp in reference_lamps:
                log.info(reference_lamp)
                if os.path.exists(reference_lamp) and os.path.isfile(reference_lamp):
                    lamp_ccd = read_fits(full_path=reference_lamp,
                                         technique=technique)
                    lamp_ccd = add_trace_info(ccd=lamp_ccd,
                                              trace_info=trace_info)

                    extracted_lamp = extraction(ccd=lamp_ccd,
                                                target_trace=trace,
                                                spatial_profile=profile,
                                                extraction_name=extraction_type)

                    new_lamp_name = "e{}".format(
                        re.sub('.fits',
                               '_{}.fits'.format(int(trace.c0.value)),
                               os.path.basename(reference_lamp)))

                    new_lamp_full_path = os.path.join(reduced_data_path, new_lamp_name)
                    write_fits(ccd=extracted_lamp,
                               full_path=new_lamp_full_path,
                               combined=False,
                               parent_file=reference_lamp,
                               overwrite=True)

                    this_target['reference_lamps'].append(new_lamp_full_path)


                else:
                    response_payload['error'] = 'Reference Lamp does not exist or is ' \
                                        'not a file'
            all_targets.append(this_target)
        response_payload['targets'] = all_targets

        return Response(response_payload)


class ApiSpectrumCalibrateView(APIView):

    def post(self, request):
        science_file = request.data.get('science_file', '')
        reference_lamps = request.data.get('reference_lamps', list())
        reduced_data_path = request.data.get('reduced_data_path')
        correlation_tolerance = request.data.get('correlation_tolerance', 15)
        technique = 'Spectroscopy'

        response_payload = {
            'api_version': __version__,
            'pipeline_version': goodman_pipeline.__version__,
            'request_data': {
                'science_file': science_file,
                'reference_lamps': reference_lamps,
                'correlation_tolerance': correlation_tolerance},
            'technique': technique,
            'error': ''
        }
        all_files = copy.deepcopy(reference_lamps)
        all_files.append(science_file)
        for _file in all_files:
            if not os.path.isfile(_file) or not os.path.exists(_file):
                response_payload['error'] = "File {} does not exist".format(science_file)
                log.error(response_payload['error'])
                return Response(response_payload,
                                status=status.HTTP_404_NOT_FOUND)

        if not os.path.exists(reduced_data_path) or not os.path.isdir(
                reduced_data_path):
            response_payload['error'] = 'reduced_data_path: {} does not exist or is '\
                                        'not a directory'.format(reduced_data_path)
            return Response(response_payload, status=status.HTTP_400_BAD_REQUEST)

        wavelength_calibration = WavelengthCalibration()


        pipeline_path = os.path.dirname(goodman_pipeline.__file__)
        reference_data_path = os.path.join(pipeline_path, 'data/ref_comp')

        ccd = read_fits(full_path=science_file, technique=technique)
        lamps = [read_fits(full_path=reference_lamp, technique=technique) for reference_lamp in reference_lamps]
        try:
            wavelength_solution = wavelength_calibration(
                ccd=ccd,
                comp_list=lamps,
                save_data_to=reduced_data_path,
                reference_data=reference_data_path,
                object_number=1,
                corr_tolerance=float(correlation_tolerance),
                plot_results=False,
                json_output=True)
            response_payload['solution'] = wavelength_solution
            return Response(response_payload)

        except NoMatchFound as error:
            response_payload['error'] = str(error)
            return Response(response_payload)
