'''
Script de prueba para la extraccion de datos del WRF para
el recurso eolico
'''
from wrfout import *
PATH = '/home/jose/Documents/wrf_data/2024_04_15_t18z/'
qii = ewrf(PATH)
lat = [-17.6290737, -17.631718, -17.62861]
lon = [-65.2842807, -65.277218, -65.28444]
qii.pm(lat, lon)
qii.forma_variables()