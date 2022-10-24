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

sys.path.insert(1,'/Users/cag6397/Desktop/EMPIC_py')
from functions_PIC import *
from functions_FDTD import *
from output import *
from definition import *
from setup import *
from parameters import *

# Abrir archivos


#######IMPORTANTE###########
#buscar una forma mas eficiente de guardar los datos de energia.
#nota: no se puede acceder a los datos de energia sino al final de la simulacion. Es importante imprimir la enegia ya que de esta forma se supervisa la evoucion del sistema.

df_energy= open('energy.txt','w')


"""## Ciclo PIC"""

i = 0 #Primera iteración
# Primera iteración ELECTRONES
v_new[:,0],v_new[:,1],v_new[:,2] = v_act[:,0],v_act[:,1],v_act[:,2]
output_particle(x_act, v_act,i, 'electrones')
v_act = velocity(x_act,v_new,ex,ey,ez,bx,by,bz,-0.5*dt,q_e,m_e, Np_e)

# Primera iteración IONES
v_new_i[:,0],v_new_i[:,1],v_new_i[:,2] = v_act_i[:,0],v_act_i[:,1],v_act_i[:,2]
output_particle(x_act_i, v_act_i, i, 'iones')
v_act_i = velocity(x_act_i,v_new_i,ex,ey,ez,bx,by,bz,-0.5*dt,q_i,m_i, Np_i)

# Guarda y avanza campos primera iteración
output_fields(ex,ey,ez,bx,by,bz,jx,jy,jz,i,divB)
bx,by,bz = FDTD_halftime(ex,ey,ez,bx,by,bz,-0.5*dt)

# Resto de iteraciones
for i in range(1, n_iteraciones):
    #PRUEBAS FDTD
    #Atmp = 0.0#np.exp(-0.5*((t0-i)/spread)**2)  # Pulso Gaussiano
    #ez[int(Nx/2),int(Ny/2)] = Atmp # Pulso Gaussiano
    #ex[int(Nx/2),int(Ny/2)]=np.exp(-i**2/20) #Superluminal waves from a source that turns on immediately.
    # ez[int(Nx/2),int(Ny/2)]=np.sin(i/5)
    # GUARDAR ENERGÍAS
    if(i % n_check == 0): # Cada n_check iteraciones guarda energias
        print('se guarda energias en iteracion', i*dt)
        energia_cinetica_electrones = energy(v_act)
        energia_cinetica_iones      = energy(v_act_i)
        energia_b = energy_fields(bx,by,bz)
        energia_e = energy_fields(ex,ey,ez)
        df_energy.write(str(i))
        df_energy.write('\t')
        df_energy.write(str(energia_cinetica_electrones))
        df_energy.write('\t')
        df_energy.write(str(energia_cinetica_iones))
        df_energy.write('\t')
        df_energy.write(str(energia_b))
        df_energy.write('\t')
        df_energy.write(str(energia_e))
        df_energy.write('\n')
        
    #incluir la energia de los campos EM

    # CICLO PRINCIPAL
    # Divergencia
    divB=divBt(bx,by,bz,divB,i)
    #guardar energia cinetica promedio del ensamble de part
    #advance E, B. E advance like: E^{n} to E^{n+1}, and B advance like: B^{n-1/2} to B^{n+1/2}
    ex,ey,ez,bx,by,bz, ex_copy, ey_copy, ez_copy, bx_copy, by_copy, bz_copy = FDTD(ex, ey, ez, bx, by, bz, jx,jy,jz,ex_copy, ey_copy, ez_copy, bx_copy, by_copy, bz_copy) #convierte los campos de la malla a DF
    #advance velocity v^{n-1/2} to v^{n+1/2}, o sea, un paso temporal
    ex_copy, ey_copy, ez_copy, bx_copy, by_copy, bz_copy = interc_to_nodes(ex, ey, ez, bx, by, bz, ex_copy, ey_copy, ez_copy, bx_copy, by_copy, bz_copy)
    v_new = velocity(x_act,v_act,ex_copy, ey_copy, ez_copy, bx_copy, by_copy, bz_copy, dt,q_e,m_e, Np_e)
    #VELOCITY IONES
    v_new_i = velocity(x_act_i,v_act_i,ex_copy, ey_copy, ez_copy, bx_copy, by_copy, bz_copy, dt, q_i, m_i, Np_i)
    #advance ELECTRONES position x^{n} to x^{n+1}, o sea, un paso temporal
    x_new, x_BC = position(x_act,v_new,x_new,x_BC, Np_e) # leapfrog
    #LEAPFROG IONES
    x_new_i, x_BC_i = position(x_act_i, v_new_i, x_new_i, x_BC_i, Np_i) # leapfrog
    # Densidad de corriente ELECTRONES
    jx = np.zeros([Nx+1, Ny+1])
    jy = np.zeros([Nx+1, Ny+1])
    jz = np.zeros([Nx+1, Ny+1])
    jx, jy = zig_zag(x_act, x_new, x_BC, jx, jy, partweight, q_e)
    # Aporte Densidad de corriente IONES -> Sobrescribe la de los electrones
    jx, jy = zig_zag(x_act_i, x_new_i, x_BC_i, jx, jy, partweight, q_i)
    # Asignar ELECTRONES antes de empezar el nuevo ciclo
    x_act[:,0],x_act[:,1],x_act[:,2] = x_BC[:,0],x_BC[:,1],x_BC[:,2]
    v_act[:,0],v_act[:,1],v_act[:,2] = v_new[:,0],v_new[:,1],v_new[:,2]
    # Asignar IONES antes de empezar el nuevo ciclo
    x_act_i[:,0], x_act_i[:,1], x_act_i[:,2] = x_BC_i[:,0], x_BC_i[:,1], x_BC_i[:,2]
    v_act_i[:,0], v_act_i[:,1], v_act_i[:,2] = v_new_i[:,0], v_new_i[:,1], v_new_i[:,2]

    # GUARDADO DE PARTÍCULAS Y CAMPOS
    if(i%n_out == 0): #Cada n_out iterations se guardan particle y campos
      #atrasa la velocidad para guardarla (la devuleve hasta tiempo {n})
      ex_copy, ey_copy, ez_copy, bx_copy, by_copy, bz_copy = interc_to_nodes(ex, ey, ez, bx, by, bz, ex_copy, ey_copy, ez_copy, bx_copy, by_copy, bz_copy)
      tmp_v_new_electrons = velocity(x_act, v_act, ex_copy, ey_copy, ez_copy, bx_copy, by_copy, bz_copy,-0.5*dt, q_e, m_e, Np_e)
      tmp_v_new_ions      = velocity(x_act_i, v_act_i, ex_copy, ey_copy, ez_copy, bx_copy, by_copy, bz_copy,-0.5*dt, q_i, m_i, Np_i)
      output_particle(x_act, tmp_v_new_electrons, i, 'electrones') # guarda electrones en tiempo {n}
      output_particle(x_act, tmp_v_new_ions, i, 'iones') # guarda IONES en tiempo {n}
      output_fields(ex_copy, ey_copy, ez_copy, bx_copy, by_copy, bz_copy,jx,jy,jz,i,divB)
      ex_tmp, ey_tmp, ez_tmp, bx_tmp, by_tmp, bz_tmp = interc_to_nodes(ex_copy, ey_copy, ez_copy, bx_copy,by_copy,bz_copy,ex_tmp, ey_tmp, ez_tmp, bx_tmp, by_tmp, bz_tmp)
      output_fields(ex_tmp,ey_tmp,ez_tmp,bx_tmp,by_tmp,bz_tmp,jx,jy,jz,i,divB,'fields_nodes')

df_energy.close()

