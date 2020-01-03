from django.urls import path

from . import views

app_name = 'api'

urlpatterns = [
    path('', views.ApiView.as_view(), name='apiview'),

   
    # path('calibrations/', views.FilesView.as_view(), name='files'),

    path('calibrations/bias/',
         views.ApiCalibrationsBiasView.as_view(),
         name='calibration_bias'),

    path('calibrations/flats/',
         views.ApiCalibrationsFlatsView.as_view(),
         name='calibration_flats'),
    path('image/overscan/',
         views.ApiImageOverscanView.as_view(),
         name='image_overscan'),
    path('image/trim/',
         views.ApiImageTrimView.as_view(),
         name='image_trim'),

    path('image/bias/',
         views.ApiImageBiasView.as_view(),
         name='image_bias'),

    path('image/flat/',
         views.ApiImageFlatView.as_view(),
         name='image_flat'),

    path('image/creject/',
         views.ApiImageCrejectView.as_view(),
         name='image_creject'),

    path('image/reduce/',
         views.ApiImageReduceView.as_view(),
         name='image_reduce'),

    path('spectrum/reduce/',
         views.ApiSpectrumReduceView.as_view(),
         name='image_reduce'),

    path('spectrum/calibrate/',
         views.ApiSpectrumCalibrateView.as_view(),
         name='image_calibrate'),
    # path('files/visualize/', views.FileVisualize.as_view(), name='fileviz'),
    # path('status/', views.ApiStatus.as_view(), name='status'),
    # path('dev/updatedb/', views.UpdateDatabase.as_view(), name='updatedb'),
]