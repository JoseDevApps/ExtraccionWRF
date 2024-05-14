'''
    Funciones para la extraccion de Salidas del WRF
'''
import numpy as np
import pandas as pd
import netCDF4
import wrf
from wrf import interplevel
import time
from glob import glob
import utm
from rich.progress import track

def ls(expr='*.*'):
    return glob(expr)

def alg_pm(x,y):
    '''
    funcion de punto medio geometrico
    '''
    n = len(x)
    if len(x)==len(y):
        sumx = np.sum(x)
        sumy = np.sum(y)
    else:
        print('Ambas listados de coordenados deben tener el mismo numero de elementos')
    return (sumx/n, sumy/n)

class ewrf():
    def __init__(self, directorio,PathOut = './',dominio=1) -> None:
        '''
        Inicializacion de la clase ewrf
        1. Especificar el directorio de salidas del WRF mediante la variable "directorio"
        2. Especificar el directorio de salida del preprocesamiento del wrfout 
            en la variable "PathOut", por defecto esta el directorio actual
        3. Especificar el dominio de las salidas a procesar mediante 
            la variable "dominio", por defecto esta en 1.

        '''
        # Declarando parametros de base
        self.directorio,self.PathOut = directorio, PathOut
        # Listado de todos los archivos del dominio declarado
        self.archivos = ls(directorio+'/*wrfout_d0'+str(dominio)+'*')
        # lectura del archivo de prueba
        self.archivo_nc = netCDF4.Dataset(self.archivos[0]) 
        print(self.archivo_nc.variables.keys())
        print('Numero de variables:',len(self.archivo_nc.variables.keys()))
    def pm(self,latitudes, longitudes):
        '''
        Punto medio geometrico de una lista de latitudes y longitudes
        '''
        self.latitudes, self.longitudes = latitudes, longitudes
        utm_zone= utm.latlon_to_zone_number(latitudes[0], longitudes[0])
        utm_zones = [utm.latlon_to_zone_number(x,y) for x,y in zip(latitudes,longitudes)]
        # Transformacion a coordenadas UTM
        utm_conv = [utm.from_latlon(x,y,force_zone_number=z) for x,y,z in zip(latitudes,longitudes,utm_zones)]
        north_cond = [x>=0 for x in latitudes]
        # Aplicacion del algoritmo
        lats = [x[0] for x in utm_conv]
        longs = [y[1] for y in utm_conv]
        self.res_UTM = alg_pm(lats, longs)
        # Convirtiendo a Lat, lon
        self.res_lat_lon = utm.to_latlon(self.res_UTM[0], self.res_UTM[1], zone_number=utm_zones[0],northern=north_cond[0])
        print([self.res_UTM, self.res_lat_lon] )
        return [self.res_UTM, self.res_lat_lon] 
    def forma_variables(self,):
        coor = wrf.ll_to_xy(self.archivo_nc, self.res_lat_lon[0], self.res_lat_lon[1], timeidx=0,
                    meta=False, as_int=True)
        #   Extrayendo velocidad de viento y direccion
        W10 = wrf.getvar(self.archivo_nc, 'uvmet10_wspd_wdir')
        #   Extrayendo Temperatura a 2m de elevacion del suelo
        T2 = wrf.getvar(self.archivo_nc, 'T2', meta=True)
        #   Extrayendo Temperatura a 2m de elevacion del suelo
        temp = wrf.getvar(self.archivo_nc, 'tc', meta=True)
        #   Extrayendo los niveles de elevacion
        z = wrf.getvar(self.archivo_nc, 'z', meta=False)
        #   Extrayendo Cloud Top Temperature
        ctt = wrf.getvar(self.archivo_nc, 'ctt')
        #   Extrayendo cloud fraction
        cf = wrf.getvar(self.archivo_nc, 'cloudfrac')
        #   Extrayendo reflectividad
        dbz = wrf.getvar(self.archivo_nc, 'dbz')
        #   Extrauendo Humedad relativa a 2m de elevacion
        rh2 = wrf.getvar(self.archivo_nc, 'rh2')
        #   Extrayendo Sea Level Pressure
        slp = wrf.getvar(self.archivo_nc, 'slp')

        # Dimensiones de las variables a ser extraidas
        Tiempos = len(self.archivos)
        print(Tiempos)
        print(W10.shape)
        print(z.shape)
        #   Referencias de forma de un archivo netcdf4
        result_shape = (Tiempos, W10.shape[0],z.shape[0], z.shape[1], z.shape[2])
        Zheight_shape = (Tiempos,z.shape[0], z.shape[1], z.shape[2])
        WindSpeed_shape = (Tiempos, W10.shape[0], W10.shape[1], W10.shape[2])
        temperature_shape = (Tiempos, W10.shape[1], W10.shape[2])
        cloud_frac_shape = (Tiempos, cf.shape[0], cf.shape[1], cf.shape[2])
        #   Reservando espacios para la extraccion del WRF
        p_final = np.empty(temperature_shape, np.float32)
        rh_final = np.empty(temperature_shape, np.float32)
        t_final = np.empty(temperature_shape, np.float32)
        W10 = np.empty(WindSpeed_shape, np.float32)
        Tiempo = ["" for x in range(0, Tiempos)]
        CTT = np.empty(temperature_shape, np.float32)
        CF = np.empty(cloud_frac_shape, np.float32)
        WWSST = np.empty(result_shape, np.float32)
        UVmet = np.empty(result_shape, np.float32)
        zz = np.empty(Zheight_shape, np.float32)
        AGLI = np.empty(Zheight_shape, np.float32)
        # Definiendo cuantas fechas existe por archivo .nc
        tiempos_p_archivo = 1
        start = time.time() # controlando los tiempos de ejecucion
        for timeidx in track(range(result_shape[0])):
            fileidx = timeidx//tiempos_p_archivo
               # Compute the file index and the time index inside the file
            file_timeidx = timeidx % tiempos_p_archivo
            f = netCDF4.Dataset(self.archivos[fileidx])
            slpw = wrf.getvar(f, 'slp', file_timeidx, meta=False)
            rh = wrf.getvar(f, 'rh2', file_timeidx, meta=False)
            t = wrf.getvar(f, 'T2', file_timeidx, meta=False)
            WS1 = wrf.getvar(f, 'uvmet10_wspd_wdir', file_timeidx, meta=False)
            Tp = wrf.getvar(f, 'times', file_timeidx, meta=False)
            WindDA = wrf.getvar(f, 'uvmet_wspd_wdir', file_timeidx, meta=False)
            UV = wrf.getvar(f, 'uvmet', file_timeidx, meta=False)  
            zh = wrf.getvar(f, "z", file_timeidx, meta=False)
            agl = wrf.getvar(f, 'height_agl', meta = False)

            WWSST[timeidx, :] = WindDA[:]
            UVmet[timeidx, :] = UV[:]
            zz[timeidx, :] = zh[:]
            AGLI[timeidx, :] =  agl[:]
            W10[timeidx, :] = WS1[:]
            p_final[timeidx, :] = slpw[:]
            rh_final[timeidx, :] = rh[:]
            t_final[timeidx, :] = t[:]
            Tiempo[timeidx] = Tp
            CTT[timeidx]  = ctt[:]
            f.close()
        end = time.time()
        print("Tiempo de ejecucion : ",end-start)
        WindDirection = np.reshape(WWSST[:,1,:,:,:], (Tiempos,z.shape[0], z.shape[1], z.shape[2]))
        WindMagnitud = np.reshape(WWSST[:,0,:,:,:], (Tiempos,z.shape[0], z.shape[1], z.shape[2]))
        WindUmet = np.reshape(UVmet[:,0,:,:,:], (Tiempos,z.shape[0], z.shape[1], z.shape[2]))
        WindVmet = np.reshape(UVmet[:,1,:,:,:], (Tiempos,z.shape[0], z.shape[1], z.shape[2]))
        Levels2 = [78, 74, 76, 57, 38]

        Vi_U = interplevel(WindUmet, AGLI, Levels2, meta=False)
        Vi_V = interplevel(WindVmet, AGLI, Levels2, meta=False)

        a_2 = Vi_U[:,:,:,:]*Vi_U[:,:,:,:]
        b_2 = Vi_V[:,:,:,:]*Vi_V[:,:,:,:]
        r = a_2 + b_2
        R = np.sqrt(r)
        # Calculo Angulo Wind Direction 78 m
        #**********************************************************
        # Angulo en coordenadas matematicas
        WD78 = np.arctan2(Vi_V[:,0,:,:], Vi_U[:,0,:,:]) * 180 / np.pi
        print(WD78.shape)
        WD78 = WD78.reshape(-1)
        print(WD78)
        # print(direccion.shape)
        for ii,tetha_i in track(enumerate(WD78), description='Convirtiendo a coordendas math...',total=100):
            if tetha_i<0.0:
                WD78[ii] = WD78[ii] + 360                 # Positive degrees
        # Conviertiendo a coordenadas terrestres donde N = 0 o 360
        # E = 90, S = 180, W = 270
        WD78 = 270 - WD78
        for ii,tetha_i in track(enumerate(WD78),description='Convirtiendo a coordendas terrestres...',total=100):
            if tetha_i<0.0:
                WD78[ii] = WD78[ii] + 360                 # Positive degrees
        np.set_printoptions(suppress=True)
        WD78 = WD78.reshape(Tiempos, z.shape[1], z.shape[2])

        # Calculo Angulo Wind Direction 74 m
        #**********************************************************
        # Angulo en coordenadas matematicas
        WD74 = np.arctan2(Vi_V[:,1,:,:], Vi_U[:,1,:,:]) * 180 / np.pi

        WD74 = WD74.reshape(-1)
        # print(direccion.shape)
        for ii,tetha_i in track(enumerate(WD74), description='Convirtiendo a coordendas math...',total=100):
            if tetha_i<0.0:
                WD74[ii] = WD74[ii] + 360                 # Positive degrees
        # Conviertiendo a coordenadas terrestres donde N = 0 o 360
        # E = 90, S = 180, W = 270
        WD74 = 270 - WD74
        for ii,tetha_i in track(enumerate(WD74), description='Convirtiendo a coordendas terrestres...',total=100):
            if tetha_i<0.0:
                WD74[ii] = WD74[ii] + 360                 # Positive degrees
        np.set_printoptions(suppress=True)
        WD74 = WD74.reshape(Tiempos, z.shape[1], z.shape[2])

        # Calculo Angulo Wind Direction 76 m
        #**********************************************************
        # Angulo en coordenadas matematicas
        WD76 = np.arctan2(Vi_V[:,2,:,:], Vi_U[:,2,:,:]) * 180 / np.pi
        WD76 = WD76.reshape(-1)
        # print(direccion.shape)
        for ii,tetha_i in track(enumerate(WD76), description='Convirtiendo a coordendas math...',total=100):
            if tetha_i<0.0:
                WD76[ii] = WD76[ii] + 360                 # Positive degrees
        # Conviertiendo a coordenadas terrestres donde N = 0 o 360
        # E = 90, S = 180, W = 270
        WD76 = 270 - WD76
        for ii,tetha_i in track(enumerate(WD76), description='Convirtiendo a coordendas terrestres...',total=100):
            if tetha_i<0.0:
                WD76[ii] = WD76[ii] + 360                 # Positive degrees
        np.set_printoptions(suppress=True)
        WD76 = WD76.reshape(Tiempos, z.shape[1], z.shape[2])


        # Calculo Angulo Wind Direction 57 m
        #**********************************************************
        # Angulo en coordenadas matematicas
        WD57 = np.arctan2(Vi_V[:,3,:,:], Vi_U[:,3,:,:]) * 180 / np.pi
        WD57 = WD57.reshape(-1)
        # print(direccion.shape)
        for ii,tetha_i in track(enumerate(WD57), description='Convirtiendo a coordendas math...',total=100):
            if tetha_i<0.0:
                WD57[ii] = WD57[ii] + 360                 # Positive degrees
        # Conviertiendo a coordenadas terrestres donde N = 0 o 360
        # E = 90, S = 180, W = 270
        WD57 = 270 - WD57
        for ii,tetha_i in track(enumerate(WD57), description='Convirtiendo a coordendas terrestres...',total=100):
            if tetha_i<0.0:
                WD57[ii] = WD57[ii] + 360                 # Positive degrees
        np.set_printoptions(suppress=True)
        WD57 = WD57.reshape(Tiempos, z.shape[1], z.shape[2])

        # Calculo Angulo Wind Direction 38 m
        #**********************************************************
        # Angulo en coordenadas matematicas
        WD38 = np.arctan2(Vi_V[:,4,:,:], Vi_U[:,4,:,:]) * 180 / np.pi
        WD38 = WD38.reshape(-1)
        # print(direccion.shape)
        for ii,tetha_i in track(enumerate(WD38), description='Convirtiendo a coordendas math...',total=100):
            if tetha_i<0.0:
                WD38[ii] = WD38[ii] + 360                 # Positive degrees
        # Conviertiendo a coordenadas terrestres donde N = 0 o 360
        # E = 90, S = 180, W = 270
        WD38 = 270 - WD38
        for ii,tetha_i in track(enumerate(WD38), description='Convirtiendo a coordendas terrestres...',total=100):
            if tetha_i<0.0:
                WD38[ii] = WD38[ii] + 360                 # Positive degrees
        np.set_printoptions(suppress=True)
        WD38 = WD38.reshape(Tiempos, z.shape[1], z.shape[2])
        # Qollpana II
        x , y = coor 
        df_QII = pd.DataFrame(
                    {'Time': Tiempo,
                    'W10s_1s': W10[:,0,y,x],
                    'Wds_1s': W10[:,1,y,x],
                    'T2s_1s': t_final[:,y,x],
                    'Rh2_1s': rh_final[:,y,x],
                    'SLP_1s': p_final[:,y,x],
                    'WS_24'+'m': WindMagnitud[:,0,y,x],
                    'WD_24'+'m': WindDirection[:,0,y,x],
                    'WS_80'+'m': WindMagnitud[:,1,y,x],
                    'WD_80'+'m': WindDirection[:,1,y,x],
                    'WS_78m': R[:,0,y,x],
                    'WS_74m': R[:,1,y,x],
                    'WS_76m': R[:,2,y,x],
                    'WS_57m': R[:,3,y,x],
                    'WS_38m': R[:,4,y,x],
                    'WD_78m': WD78[:,y,x],
                    'WD_74m': WD74[:,y,x],
                    'WD_76m': WD76[:,y,x],
                    'WD_57m': WD57[:,y,x],
                    'WD_38m': WD38[:,y,x]
                    },index=None
        )
        df_QII = df_QII[['Time', 'W10s_1s', 'Wds_1s','T2s_1s','Rh2_1s','SLP_1s','WS_24m','WD_24m','WS_80m','WD_80m','WS_78m', 'WS_74m', 'WS_76m', 'WS_57m', 'WS_38m','WD_78m', 'WD_74m', 'WD_76m', 'WD_57m', 'WD_38m']]          # set order of columns
        print(df_QII)
        # Convertir a UTC-4
        df_QII['Time'] = df_QII['Time'] - pd.Timedelta(hours=4)
        df_QII = df_QII.sort_values(by='Time')
        df_QII.to_csv(self.PathOut +'QII'+ '.csv',index=False)