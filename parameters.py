# coding=utf-8
import math
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from numpy import random
import pandas as pd
import random
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

# paso temporal y espacial
dx            = 1.0
dy            = 1.0
dt            = 0.1
Nx            = 64
Ny            = 64
Nz            = 1
Xmin          = 0
Xmax          = Nx
Ymin          = 0
Ymax          = Ny
Lx            = Xmax
Ly            = Ymax

# Parámetros del ciclo (iteraciones)
t_final       = 200
n_check       = 1.0 # Cada cuanto se guarda informacion de energia
n_out         = 5.0 # Cada cuanto se guarda outputs (campos y particulas)

#####IMPORTANTE#################
#crear variables para campos y particular. Nota: flexibilidad en guardar archivos en distintos timpos, es muuy pesado guardar la informacion de particulas con mucha frequencia.
#n_out_particles
#n_out_fields

n_iteraciones = int(t_final/dt)
n_check = int(n_check/dt)
n_out   = int(n_out/dt)


# Parámetros de las partículas en la malla
Npart_cell    = 10 # Número inicial de PARTÍCULAS en cada celda
Np_i          = Npart_cell*Nx*Ny # Número total de electrones
Np_e          = Npart_cell*Nx*Ny # Número total de iones
Np            = Np_e # Número total de partículas
partweight    = 1/Npart_cell # Peso de cada partícula por celda

# PARCMETROS DEL PLASMA (Gaussian cgs units)
m_e           = 1.0 #masa del electrC3n
c             = 1.0 #velocidad de la luz
q_0           = 1.0 #carga elC)ctrica de referencia
kb            = 1.0 #Constante de Boltzmann
n0            = 1.0 #Densidad de electrones
eps           = 1.0 #permitividad del vacio (Faraday/m)
mi_me         = 100 #razC3n de masa ion/ele
m_i           = mi_me*m_e #masa iones
q_e           = -1*q_0 #carga del electrC3n
q_i           = 1*q_0 #carga iC3n (protC3n)
courant_num   = (c*dt)/dx
beta_e        = 8.0 # magnetic electron rest enegy ratio
beta_i        = 2.0 # magnetic ion rest enegy ratio

#Parametros del plasma
w_pi          = np.sqrt((n0*q_0**2)/m_i) #Frecuencia de los iones en el plasma
w_pe          = np.sqrt((n0*q_0**2)/m_e) #Frecuencia de los electrones en el plasma
B_0           = m_e*w_pe/q_0 # PerturbaciC3n inicial Campo magnC)tico
w_ce          = q_0*B_0/m_e*c #electron gyrofrequency
w_ci          = q_0*B_0/m_i*c #ion gyrofrequency
Ti            = (B_0**2 * beta_i)/(2*n0*kb) #Temperatura de los iones (K)
Te            = (B_0**2 * beta_e)/(2*n0*kb) #Temperatura de los electrones (K)
vth_e         = (kb*Te/m_e)**0.5 # Velocidad térmica electrones
vth_i         = (kb*Ti/m_i)**0.5 # Velocidad térmica iones
lambda_d      = ((kb*Te)/(n0*q_0**2))**0.5 # Longitud de Debye
r_e           = vth_e/w_ce # Giroradio del electrón
r_i           = vth_i/w_ci # Giroradio del ion
d_e           = c/w_pe # distancia inercial del electrón
d_i           = c/w_pi # distancia inercial del protón

print('Parametros del plasma (no-normalizados )')
print('escalas espaciales: ')
print('l_d=',lambda_d, 'r_e=',r_e,'d_e',d_e)
print('r_i=',r_i,'d_i=',d_i)
print('escalas temporales: ')
print('w_pe=',w_pe,'w_ce=',w_ce)
print('w_pi=',w_pi,'w_ci=',w_ci)
print('velocidades:')
print('vth_e=',vth_e,'vth_i=',vth_i)

# Parámetros de renomalización
w_0           = w_pe
x_0           = c/w_pe
v_0           = c
E0            = m_e*c*w_0/q_0
N0            = (eps*m_e*w_0*w_0)/(q_0*q_0)
J0            = c*q_0*N0

# Renormalización
wci          = w_ci/w_0
wce          = w_ce/w_0
wpi          = w_pi/w_0
wpe          = w_pe/w_0
vthi         = vth_i/v_0
vthe         = vth_e/v_0
di           = d_i/x_0
de           = d_e/x_0
ri           = r_i/x_0
re           = r_e/x_0
lambdad      = lambda_d/x_0
print ('\n')
print('Parametros del plasma (normalizados )')
print('escalas espaciales: ')
print('l_d=',lambdad, 'r_e=',re,'r_i=',ri)
print('escalas temporales: ')
print('w_pe=',wpe,'w_pi=',wpi)
print('velocidades:')
print('vth_e=',vthe,'vth_i=',vthi)#imprimir escalas temporales
print ('\n\n\n')

# Pendiente
dim           = 3

# Para la distribución gaussiana, qué poner para que quede solamente la vth
mu            = 1
sigma         = 1
Vth           = mu

