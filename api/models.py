import os

from django.db import models
from django.contrib.auth.models import User


# class ApiSettings(models.Model):
#     user = models.ForeignKey(User, on_delete=models.CASCADE)
#     data_path = models.CharField(max_length=200)
#     file_pattern = models.CharField(max_length=20)
#     sub_folder = models.CharField(max_length=200)


# class ApiSessions(models.Model):
#     created = models.DateTimeField(auto_now_add=True)
#     observer_name = models.CharField(max_length=50)
#     # observing_technique
#     raw_data_path = models.CharField(max_length=200)
#     reduced_data_path = models.CharField(max_length=200)
#     overscan_region = models.CharField(max_length=50)
#     trim_section = models.CharField(max_length=50)
#     slit_trim_section = models.CharField(max_length=50)
#     master_bias = models.CharField(max_length=200)
#     master_flat = models.CharField(max_length=200)


# class Calibrations(models.Model):
#
#     session = models.ForeignKey(ApiSessions, on_delete=models.CASCADE)
#     file_name = models.CharField(max_length=200, unique=True)
#     full_path = models.CharField(max_length=200)
#     obstype = models.CharField(max_length=10)
#     gain = models.CharField(max_length=10)
#     read_noise = models.CharField(max_length=10)
#     filter = models.CharField(max_length=200, blank=True, null=True)
#     filter_2 = models.CharField(max_length=200, blank=True, null=True)
#     grating = models.CharField(max_length=200, blank=True, null=True)
#     slit = models.CharField(max_length=200, blank=True, null=True)
#     cam_targ = models.CharField(max_length=200, blank=True, null=True)
#     grt_targ = models.CharField(max_length=200, blank=True, null=True)
#     roi = models.CharField(max_length=200, blank=True, null=True)
#
#     def as_json(self):
#         return {"session": self.session.id,
#                 "file_name": self.file_name,
#                 "full_path": self.full_path,
#                 "obstype": self.obstype,
#                 "gain": self.gain,
#                 "read_noise": self.read_noise,
#                 "filter": self.filter,
#                 "filter2": self.filter_2,
#                 "grating": self.grating,
#                 "slit": self.slit,
#                 "cam_targ": self.cam_targ,
#                 "grt_targ": self.grt_targ,
#                 "roi": self.roi
#                 }
#
#     def file_exist(self):
#         return os.path.exists(os.path.join(self.full_path, self.file_name))

#
# class ObservedFile(models.Model):
#     file_name = models.CharField(max_length=200, unique=True)
#     file_status = models.CharField(max_length=20, default='EXIST')
#     object = models.CharField(max_length=200)
#     obs_date = models.DateField()
#     obs_time = models.DateTimeField()
