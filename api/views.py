from asgiref.sync import async_to_sync
from astropy.io import fits
from ccdproc import CCDData
import ccdproc
import matplotlib.pyplot as plt
import glob
import json
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

setup_logging()


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
        raise FileNotFoundError


def find_matching_calibration(full_path, calibration_type, search_path):
    if os.path.isdir(os.path.dirname(full_path)) and os.path.exists(full_path):
        if not os.path.isabs(full_path):
            full_path = os.path.abspath(full_path)

        header = fits.getheader(full_path)
        ifc = ImageFileCollection(location=os.path.abspath(search_path),
                                  keywords=KEYWORDS)
        pd_ifc = ifc.summary.to_pandas()

        if calibration_type == 'BIAS':
            filtered = pd_ifc[
                ((pd_ifc['OBSTYPE'] == calibration_type) &
                 (pd_ifc['GAIN'] == header['GAIN']) &
                 (pd_ifc['RDNOISE'] == header['RDNOISE']) &
                 (pd_ifc['ROI'] == header['ROI'])
                 )]
            bias_list = filtered.file.tolist()

            if len(bias_list) == 1:
                return bias_list[0]
            elif len(bias_list) > 1:
                return bias_list[-1]
            else:
                raise FileNotFoundError

        elif calibration_type == 'FLAT':
            filtered = pd_ifc[
                (((pd_ifc['OBSTYPE'] == 'LAMPFLAT') |
                  (pd_ifc['OBSTYPE'] == 'FLAT')) &
                 (pd_ifc['FILTER'] == header['FILTER']) &
                 (pd_ifc['FILTER2'] == header['FILTER2']) &
                 (pd_ifc['GRATING'] == header['GRATING']) &
                 (pd_ifc['SLIT'] == header['SLIT']) &
                 (pd_ifc['CAM_TARG'] == header['CAM_TARG']) &
                 (pd_ifc['GRT_TARG'] == header['GRT_TARG']) &
                 (pd_ifc['GAIN'] == header['GAIN']) &
                 (pd_ifc['RDNOISE'] == header['RDNOISE']) &
                 (pd_ifc['ROI'] == header['ROI'])
                 )]
            flat_list = filtered.file.tolist()
            if len(flat_list) == 1:
                return flat_list[0]
            elif len(flat_list) > 1:
                return flat_list[-1]
            else:
                raise FileNotFoundError

    else:
        raise FileNotFoundError


class ApiView(APIView):
    def get(self, request):
        content = {'detail': 'Data Reduction Goodman Pipeline API'}
        return Response(content)


class ApiCalibrationsBiasView(APIView):
    
    def get(self, request):
        file_name = request.GET.get('file', '')
        raw_data_path = request.GET.get('raw_data_path')
        reduced_data_path = request.GET.get('reduced_data_path')

        if file_name != '':
            file_full_path = os.path.join(reduced_data_path, file_name)
            if os.path.exists(file_full_path):
                try:
                    master_bias_file = find_matching_calibration(
                        full_path=file_full_path,
                        calibration_type='BIAS',
                        search_path=reduced_data_path)

                    return Response({
                        'api_version': __version__,
                        'pipeline_version': goodman_pipeline.__version__,
                        "master_bias": master_bias_file})
                except FileNotFoundError:
                    return Response({
                        'api_version': __version__,
                        'pipeline_version': goodman_pipeline.__version__,
                        'error': 'No master bias was found for {}'
                                 ''.format(file_name)})

            else:
                return Response({
                    'api_version': __version__,
                    'pipeline_version': goodman_pipeline.__version__,
                    "error": "File {} not found".format(file_name)},
                    status=status.HTTP_404_NOT_FOUND)
        else:

            all_master_bias = auto_discover_files(full_path=reduced_data_path,
                                                  file_type='BIAS')
            if len(all_master_bias) == 1:
                return Response({
                    'api_version': __version__,
                    'pipeline_version': goodman_pipeline.__version__,
                    'master_bias': all_master_bias})
            else:
                return Response({
                    'api_version': __version__,
                    'pipeline_version': goodman_pipeline.__version__,
                    'error': "Unable to find master bias files"})

    def post(self, request):
        auto_discover = request.data.get('auto_discover')
        raw_data_path = request.data.get('raw_data_path')
        reduced_data_path = request.data.get('reduced_data_path')
        overscan_region = request.data.get('overscan_region')
        trim_section = request.data.get('trim_section')

        if auto_discover == 'true':
            try:
                file_list = auto_discover_files(full_path=raw_data_path,
                                                file_type='BIAS')
            except FileNotFoundError:
                response_payload = {
                    'api_version': __version__,
                    'pipeline_version': goodman_pipeline.__version__,
                    'error': 'Directory {} could not '
                             'be found'.format(raw_data_path)}
                return Response(response_payload,
                                status=status.HTTP_404_NOT_FOUND)
        else:
            file_list = request.data.get('file_list')

        if file_list != []:
            all_master_bias = []
            for sub_list in file_list:
                if isinstance(sub_list, list):
                    start_time = time.time()
                    if overscan_region == "" or trim_section == "":
                        sample_file = os.path.join(raw_data_path, sub_list[0])

                        if trim_section == "":
                            trim_section = define_trim_section(
                                sample_image=sample_file,
                                technique='Spectroscopy')

                        if overscan_region == "":
                            overscan_region = get_overscan_region(
                                sample_image=sample_file,
                                technique="Spectroscopy")
                    master_bias, master_bias_name = create_master_bias(
                        bias_files=sub_list,
                        raw_data=raw_data_path,
                        reduced_data=reduced_data_path,
                        technique='Spectroscopy',
                        overscan_region=overscan_region,
                        trim_section=trim_section)
                    elapsed = time.time() - start_time

                    bias_name = re.sub('.fits', '', os.path.basename(master_bias_name))

                    all_master_bias.append({bias_name: {
                        'file_name': os.path.basename(master_bias_name),
                        'full_path': os.path.dirname(master_bias_name),
                        'file_list': sub_list,
                        'elapsed_time': elapsed}})

            response_payload = {
                'api_version': __version__,
                'pipeline_version': goodman_pipeline.__version__,
                'master_bias': all_master_bias}

            return Response(response_payload)
        else:
            return Response({
                'api_version': __version__,
                'pipeline_version': goodman_pipeline.__version__,
                'error': "Unable to create master flats"})


class ApiCalibrationsFlatsView(APIView):

    def get(self, request):

        sample_file = request.GET.get('file', '')
        raw_data_path = request.GET.get('raw_data_path', '')
        reduced_data_path = request.GET.get('reduced_data_path', '')

        if sample_file != '':
            full_path = os.path.join(raw_data_path, sample_file)
            if os.path.exists(full_path):
                master_flat_file = find_matching_calibration(
                    full_path=full_path,
                    calibration_type='FLAT',
                    search_path=reduced_data_path)

                return Response({
                    'api_version': __version__,
                    'pipeline_version': goodman_pipeline.__version__,
                    'master_flat': master_flat_file})

            else:
                return Response({
                    'api_version': __version__,
                    'pipeline_version': goodman_pipeline.__version__,
                    'error': 'File {} does not exist'.format(sample_file)})
        else:
            all_flats = auto_discover_files(full_path=reduced_data_path,
                                            file_type='FLAT')

        return Response({'api_version': __version__,
                         'pipeline_version': goodman_pipeline.__version__,
                         'master_flat': all_flats})

    def post(self, request):
        auto_discover = request.data.get('auto_discover')
        saturation_threshold = request.data.get('saturation_threshold')
        raw_data_path = request.data.get('raw_data_path')
        reduced_data_path = request.data.get('reduced_data_path')
        overscan_region = request.data.get('overscan_region')
        trim_section = request.data.get('trim_section')
        normalize_method = request.data.get('normalize_method', 'simple')
        polynomial_order = request.data.get('polynomial_order', 15)

        response_payload = {'api_version': __version__,
                            'pipeline_version': goodman_pipeline.__version__}

        ignore_bias = request.data.get('ignore_bias')
        if ignore_bias == 'false':
            ignore_bias = False
        elif ignore_bias == 'true':
            ignore_bias = True

        if not ignore_bias:
            master_bias_name = request.data.get('master_bias')
            if not os.path.exists(os.path.join(reduced_data_path,
                                               master_bias_name)):
                response_payload['error'] = "Unable to locate master bias file "\
                                            "{}".format(master_bias_name)
                return Response(response_payload,
                                status=status.HTTP_404_NOT_FOUND)
        else:
            master_bias_name = ''

        if auto_discover == 'true':
            all_files = auto_discover_files(full_path=raw_data_path,
                                            file_type='FLAT')
        else:
            all_files = request.data.get('file_list')

        all_flats = []
        for _file_list in all_files:
            sample_file = os.path.join(raw_data_path,
                                       _file_list[0])

            if overscan_region == "":
                overscan_region = get_overscan_region(sample_image=sample_file,
                                                      technique='Spectroscopy')
                response_payload['overscan_region'] = overscan_region
            if trim_section == "":
                trim_section = define_trim_section(sample_image=sample_file,
                                                   technique='Spectroscopy')
                response_payload['trim_section'] = trim_section

            file_header = fits.getheader(sample_file)

            new_master_flat_name = "master_flat_{}_{}.fits".format(
                re.sub(' ', '_', file_header['WAVMODE']),
                re.sub('[<> ]', '', file_header['FILTER2']))

            try:
                master_flat, master_flat_name = create_master_flats(
                    flat_files=_file_list,
                    raw_data=raw_data_path,
                    reduced_data=reduced_data_path,
                    technique="Spectroscopy",
                    overscan_region=overscan_region,
                    trim_section=trim_section,
                    master_bias_name=master_bias_name,
                    new_master_flat_name=new_master_flat_name,
                    saturation=saturation_threshold,
                    ignore_bias=ignore_bias)

                norm_master,  norm_master_name = normalize_master_flat(
                    master=master_flat,
                    name=master_flat_name,
                    method=normalize_method,
                    order=int(polynomial_order))

                entry_name = re.sub('.fits', '',
                                    os.path.basename(master_flat_name))
                all_flats.append(
                    {entry_name: {
                        'master_flat': os.path.basename(master_flat_name),
                        'norm_master_flat': os.path.basename(norm_master_name),
                        'full_path': os.path.dirname(master_flat_name),
                        'file_list': _file_list}})
            except ValueError as error:
                response_payload['error'] = str(error)
                return Response(response_payload)

        response_payload['master_flat'] = all_flats
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
            if master_bias_file == '':
                try:
                    master_bias_file = find_matching_calibration(
                        full_path=full_path,
                        calibration_type='BIAS',
                        search_path=reduced_data_path)
                    master_bias = read_fits(
                        full_path=os.path.join(
                            reduced_data_path,
                            master_bias_file),
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
                except FileNotFoundError:
                    return Response({'api_version': __version__,
                             'pipeline_version': goodman_pipeline.__version__,
                        'error': "Unable to find a suitable master "
                                 "bias for file {}".format(parent_file)})
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

    def get(self, request):
        return Response({'api_version': __version__,
                         'pipeline_version': goodman_pipeline.__version__,
                         'error': 'not implemented'})

    def post(self, request):
        parent_file = request.data.get('file', '')
        if parent_file == '':
            return Response({'api_version': __version__,
                             'pipeline_version': goodman_pipeline.__version__,
                             'error': 'Must provide a valid file name'})

        reduced_data_path = request.data.get('reduced_data_path')
        method = request.data.get('method', 'default')

        if os.path.exists(reduced_data_path) and os.path.isdir(
                reduced_data_path):
            full_path = os.path.join(reduced_data_path, parent_file)
            if not os.path.exists(full_path) or not os.path.isfile(full_path):
                return Response(
                    {'api_version': __version__,
                     'pipeline_version': goodman_pipeline.__version__,
                     'error': 'file {} does not exist'.format(parent_file)})

            ccd = read_fits(full_path=full_path, technique='Spectroscopy')

            ccd, astroscrappy_lacosmic(
                ccd=ccd,
                red_path=reduced_data_path,
                save_mask=False)

            new_file_name = "c{}".format(parent_file)
            full_path = os.path.join(reduced_data_path, new_file_name)

            write_fits(ccd=ccd, full_path=full_path)

            return Response({'api_version': __version__,
                             'pipeline_version': goodman_pipeline.__version__,'crejected': os.path.basename(full_path)})

        else:
            return Response({'api_version': __version__,
                             'pipeline_version': goodman_pipeline.__version__,
                'error': 'folder {} does not exist'.format(
                    reduced_data_path)})


class ApiImageReduceView(APIView):

    def post(self, request):
        file_name = request.data.get('file', '')
        if file_name == '':
            return Response({'api_version': __version__,
                             'pipeline_version': goodman_pipeline.__version__,
                             'error': 'please provide a valid file name.'})
        raw_data_path = request.data.get('raw_data_path')
        reduced_data_path = request.data.get('reduced_data_path')
        technique = request.data.get('technique', 'Spectroscopy')

        if not os.path.exists(raw_data_path) or not os.path.isdir(raw_data_path):
            return Response({'api_version': __version__,
                             'pipeline_version': goodman_pipeline.__version__,
                             'error': 'raw_data_path does not exist or is '
                                      'not a directory'})
        if not os.path.exists(reduced_data_path) or not os.path.isdir(reduced_data_path):
            return Response({'api_version': __version__,
                             'pipeline_version': goodman_pipeline.__version__,
                             'error': 'reduced_data_path does not exist or is '
                                      'not a directory'})

        full_path = os.path.join(raw_data_path, file_name)

        if not os.path.exists(full_path) or not os.path.isfile(full_path):
            return Response({'api_version': __version__,
                             'pipeline_version': goodman_pipeline.__version__,
                             'error': "{} does not exist in raw_data_path or "
                                      "is not a file".format(file_name)})

        overscan_region = get_overscan_region(sample_image=full_path,
                                              technique=technique)

        trim_section = define_trim_section(sample_image=full_path,
                                           technique=technique)

        try:
            master_bias_name = find_matching_calibration(
                full_path=full_path,
                calibration_type='BIAS',
                search_path=reduced_data_path)
            master_bias = read_fits(full_path=os.path.join(reduced_data_path,
                                                           master_bias_name),
                                    technique=technique)

            master_bias = image_trim(ccd=master_bias,
                                     trim_section=trim_section,
                                     trim_type='trimsec',
                                     add_keyword=False)

        except FileNotFoundError:
            return Response({'api_version': __version__,
                             'pipeline_version': goodman_pipeline.__version__,
                             'error': 'Master bias not found for file '
                                      '{}'.format(file_name)})

        try:
            master_flat_name = find_matching_calibration(
                full_path=full_path,
                calibration_type='FLAT',
                search_path=reduced_data_path)
            master_flat = read_fits(full_path=os.path.join(reduced_data_path,
                                                           master_flat_name),
                                    technique=technique)

            master_flat = image_trim(ccd=master_flat,
                                     trim_section=trim_section,
                                     trim_type='trimsec',
                                     add_keyword=False)

        except FileNotFoundError:
            return Response({'api_version': __version__,
                             'pipeline_version': goodman_pipeline.__version__,
                             'error': 'Master flat not found for file '
                                      '{}'.format(file_name)})

        slit_trim_section = get_slit_trim_section(master_flat=master_flat)

        master_bias = image_trim(ccd=master_bias,
                                 trim_section=slit_trim_section,
                                 trim_type='slit',
                                 add_keyword=False)

        master_flat = image_trim(ccd=master_flat,
                                 trim_section=slit_trim_section,
                                 trim_type='slit',
                                 add_keyword=False)

        prefix = '_'

        ccd = read_fits(full_path=full_path, technique=technique)

        ccd = image_overscan(ccd=ccd,
                             overscan_region=overscan_region,
                             add_keyword=False)
        prefix = 'o' + prefix

        # save intermediate?

        ccd = image_trim(ccd=ccd,
                         trim_section=trim_section,
                         trim_type='trimsec')
        prefix = 't' + prefix

        ccd = image_trim(ccd=ccd,
                         trim_section=slit_trim_section,
                         trim_type='slit')

        prefix = 's' + prefix

        ccd = ccdproc.subtract_bias(ccd=ccd,
                                    master=master_bias,
                                    add_keyword=False)

        prefix = 'z' + prefix

        ccd = ccdproc.flat_correct(ccd=ccd,
                                   flat=master_flat,
                                   add_keyword=False)

        prefix = 'f' + prefix

        ccd = astroscrappy_lacosmic(ccd=ccd,
                                    red_path=reduced_data_path,
                                    save_mask=False)

        prefix = 'c' + prefix


        new_name = prefix + file_name
        new_full_path = os.path.join(reduced_data_path, new_name)

        ccd.write(new_full_path, overwrite=True)

        return Response({'api_version': __version__,
                         'pipeline_version': goodman_pipeline.__version__,
                         'overscan_region': overscan_region,
                         'trim_section': trim_section,
                         'slit_trim_section': slit_trim_section,
                         'master_bias': master_bias_name,
                         'master_flat': master_flat_name,
                         'reduced': new_name})


class ApiSpectrumReduceView(APIView):

    def post(self, request):
        response = {}
        file_name = request.data.get('file', '')
        if file_name == '':
            return Response({'api_version': __version__,
                             'pipeline_version': goodman_pipeline.__version__,
                             'error': 'please provide a valid file name.'})
        data_path = request.data.get('data_path')

        reference_lamp = request.data.get('reference_lamp')
        if reference_lamp == '':
            response['warning'] = 'Reference lamp not provided'

        fit_model = request.data.get('fit_model', 'moffat')
        background_threshold = request.data.get('background_threshold', 1)
        max_identify = request.data.get('max_identify', 3)
        extraction_type = request.data.get('extraction_type', 'fractional')
        technique = 'Spectroscopy'

        if not os.path.exists(data_path) or not os.path.isdir(
                data_path):
            return Response({'api_version': __version__,
                             'pipeline_version': goodman_pipeline.__version__,
                             'error': 'data_path: {} does not exist or is '
                                      'not a directory'.format(data_path)})

        full_path = os.path.join(data_path, file_name)

        if not os.path.exists(full_path) or not os.path.isfile(full_path):
            return Response({'api_version': __version__,
                             'pipeline_version': goodman_pipeline.__version__,
                             'error': "{} does not exist in raw_data_path or "
                                      "is not a file".format(file_name)})

        ccd = read_fits(full_path=full_path, technique=technique)

        targets = identify_targets(ccd=ccd,
                                   fit_model=fit_model,
                                   background_threshold=background_threshold,
                                   nfind=max_identify)
        if len(targets) == 0:
            return Response({'api_version': __version__,
                             'pipeline_version': goodman_pipeline.__version__,
                             'error': 'Unable to identify targets in '
                                      '{}'.format(file_name)})

        traces = trace_targets(ccd=ccd,
                               target_list=targets,
                               sampling_step=5,
                               pol_deg=2,
                               nfwhm=5,
                               plots=False)

        i = 1
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
                       file_name))

            new_full_path = os.path.join(data_path, new_file_name)
            write_fits(ccd=extracted_ccd,
                       full_path=new_full_path,
                       combined=False,
                       parent_file=file_name,
                       overwrite=True)

            response[target_name] = {'target_center': trace.c0.value,
                                     'file_name': new_file_name,
                                     'trace_rms_error': trace_info['GSP_TERR'][0]}
            lamp_full_path = os.path.join(data_path, reference_lamp)

            if os.path.exists(lamp_full_path) and os.path.isfile(lamp_full_path):
                lamp_ccd = read_fits(full_path=lamp_full_path,
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
                           reference_lamp))

                new_lamp_full_path = os.path.join(data_path, new_lamp_name)
                write_fits(ccd=extracted_lamp,
                           full_path=new_lamp_full_path,
                           combined=False,
                           parent_file=reference_lamp,
                           overwrite=True)

                response[target_name]['reference_lamp'] = new_lamp_name

            else:
                response['error'] = 'Reference Lamp does not exist or is ' \
                                    'not a file'

        return Response(response)


class ApiSpectrumCalibrateView(APIView):

    def post(self, request):
        file_name = request.data.get('file', '')
        if file_name == '':
            return Response({'api_version': __version__,
                             'pipeline_version': goodman_pipeline.__version__,
                             'error': 'please provide a valid file name.'})

        reference_lamp = request.data.get('reference_lamp')
        if reference_lamp == '':
            return Response({'api_version': __version__,
                             'pipeline_version': goodman_pipeline.__version__,
                             'error': 'Reference lamp file not provided'})

        data_path = request.data.get('data_path')
        if not os.path.exists(data_path) or not os.path.isdir(
                data_path):
            return Response({'api_version': __version__,
                             'pipeline_version': goodman_pipeline.__version__,
                             'error': 'data_path: {} does not exist or is '
                                      'not a directory'.format(data_path)})
        correlation_tolerance = request.data.get('correlation_tolerance', 15)

        file_full_path = os.path.join(data_path, file_name)
        if not os.path.exists(file_full_path) or not os.path.isfile(file_full_path):
            return Response({'api_version': __version__,
                             'pipeline_version': goodman_pipeline.__version__,
                             'error': 'File {} does not exist or is not an '
                                      'actual file'.format(file_name)})

        lamp_full_path = os.path.join(data_path, reference_lamp)
        if not os.path.exists(lamp_full_path) or not os.path.isfile(lamp_full_path):
            return Response({'api_version': __version__,
                             'pipeline_version': goodman_pipeline.__version__,
                             'error': "Reference lamp {} does not exist or is "
                                      "not a file".format(reference_lamp)})
        wavelength_calibration = WavelengthCalibration()

        pipeline_path = os.path.dirname(goodman_pipeline.__file__)
        reference_data_path = os.path.join(pipeline_path, 'data/ref_comp')

        ccd = read_fits(full_path=file_full_path, technique='Spectroscopy')
        lamp = read_fits(full_path=lamp_full_path, technique='Spectroscopy')
        try:
            wavelength_solution = wavelength_calibration(
                ccd=ccd,
                comp_list=[lamp],
                save_data_to=data_path,
                reference_data=reference_data_path,
                object_number=1,
                corr_tolerance=float(correlation_tolerance),
                plot_results=False,
                json_output=True)
            return Response({'api_version': __version__,
                             'pipeline_version': goodman_pipeline.__version__,
                             'solution': wavelength_solution})

        except NoMatchFound as error:
            return Response({'api_version': __version__,
                             'pipeline_version': goodman_pipeline.__version__,
                             'error': str(error)})

# def api_view(request):
#     return render(request=request, template_name='api/api.html')


# class FilesView(ApiView):
# 
#     permission_classes = (IsAuthenticated,)
# 
#     def get(self, request):
# 
#         auth_token = request.headers['Authorization']
# 
#         settings = ApiSettings.objects.all()[0]
#         
#         files_in_record = ObservedFile.objects.all().values_list('file_name',
#                                                                 flat=True)
# 
#         print(files_in_record)
# 
#         result = sorted(glob.glob(os.path.join(settings.data_path,
#                                                settings.file_pattern)))
# 
#         file_list = [os.path.basename(_file) for _file in result]
# 
#         new_files = [_file for _file in file_list if _file not in files_in_record]
#         missing_files = [_file for _file in files_in_record if _file not in file_list]
# 
#         if new_files != []:
#             for _file in new_files:
#                 header_response = requests.get(
#                     'http://localhost:8000/api/files/header/',
#                     params={'filename': _file},
#                     headers={'Authorization': auth_token})
#                 header_data = json.loads(header_response.text)
#                 header = header_data['header']
# 
#                 _f = ObservedFile(file_name=_file,
#                                   object=header['OBJECT'],
#                                   obs_date=header['DATE'],
#                                   obs_time=header['DATE-OBS'])
#                 _f.save()
# 
# 
#         
#         files_in_record = ObservedFile.objects.all().order_by('-id')
# 
# 
#         all_files = []
#         #
#         for _file in files_in_record:
#             res = dict({
#                 'file_name': _file.file_name,
#                 'object': _file.object,
#                 'obs_date': _file.obs_date,
#                 'obs_time': _file.obs_time
#             })
# 
#             all_files.append(res)
# 
#         print(all_files)
# 
#         response_data = dict()
# 
#         response_data['data_path'] = settings.data_path
#         response_data['file_pattern'] = settings.file_pattern
#         # response_data['query_set'] = files_in_record
#         response_data['file_list'] = file_list
#         response_data['new_files'] = new_files
#         response_data['missing_files'] = missing_files
# 
#         return Response(response_data)
# 
# 
# class FileHeader(APIView):
# 
#     permission_classes = (IsAuthenticated,)
# 
#     def get(self, request):
#         filename = request.GET.get('filename')
#         response_data = dict({'filename': '',
#                               'header': '',
#                               'error': ''})
# 
#         if filename is None:
#             response_data['error'] = "No filename querystring was provided"
#         else:
#             response_data['filename'] = filename
#             settings = ApiSettings.objects.all()[0]
# 
#             if os.path.isfile(os.path.join(settings.data_path, filename)):
#                 fits_header = fits.getheader(os.path.join(settings.data_path,
#                                                           filename))
# 
#                 header = dict()
# 
#                 for key in fits_header.keys():
#                     if key != '':
#                         header[key] = str(fits_header[key])
# 
#                 response_data['header'] = header
# 
#             else:
#                 response_data['error'] = 'File {:s} not found'.format(filename)
# 
#         return Response(response_data)
# 
# 
# class FileVisualize(APIView):
# 
#     permission_classes = (IsAuthenticated,)
# 
#     def get(self, request):
#         filename = request.GET.get('filename')
# 
#         response_data = dict({'filename': '',
#                               'error': ''})
# 
#         if filename is None:
#             response_data['error'] = "No filename querystring was provided"
#         else:
#             response_data['filename'] = filename
#             settings = ApiSettings.objects.all()[0]
#             if os.path.isfile(os.path.join(settings.data_path, filename)):
#                 response_data['filename'] = filename
# 
#                 plot_file = _make_plot(full_path=os.path.join(settings.data_path,
#                                                               filename),
#                                        static_path='live/static/plots')
#                 response_data['plot'] = plot_file
# 
#             else:
#                 response_data['error'] = 'File {:s} not found'.format(filename)
# 
#         return Response(response_data)
