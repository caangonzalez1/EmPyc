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
"""# Funciones Graficos#"""

def graficar_particula(especie, i,typecl):
#  fig, ax = plt.subplots()
#  plt.title("Trayectoria de la partícula")
#  plt.xlabel("x")
#  plt.ylabel("y")
  # Particle Trajectory
  lp = [] #np.array([])
  lv = [] #np.array([])
  ii= 50 #particular N
  print(ii)

  for i in range(0,50,5):
    file_name = especie + str(i) +'.pkl'
    print(file_name)
    df = pd.read_pickle(file_name)
    tmp1 = np.array(df.iloc[ii,[0,1,2]])
    tmp2 = np.array(df.iloc[ii,[3,4,5]])
    # print(lp)
    # input()
    lp.append(tmp1) #lp = np.append(lp, tmp1, axis=i) #
    lv.append(tmp2) #lv = np.append(lv, tmp2, axis=i)
  lp = np.array(lp)
  lv = np.array(lv)

  for i in range(len(lp)):
    plt.plot(lp[i,0], lp[i,1],typecl)
    plt.plot()

  # i = 10
  # file_name = especie + str(i) +'.pkl'
  # print(file_name)
  # df = pd.read_pickle(file_name)
  # x_act = np.array(df.iloc[:,[0,1,2]])
  # v_act = np.array(df.iloc[:,[3,4,5]])

  # v_new = velocity(x_act,v_act,ex,ey,ez,bx,by,bz,dt, q_e, m_e)
  # x_new,x_BC = position(x_act,v_new,x_new,x_BC)
  # jx,jy = zig_zag(x_act,x_new,x_BC,jx,jy,partweight)
  # jx
  # jy

  # plt.figure(13);plt.clf()

  # plt.subplot(1,2,1)
  # plt.rcParams['figure.figsize'] = [14, 6]
  # plt.rcParams['figure.dpi'] = 50 # 200 e.g. is really fine, but slower
  # plt.title("Densidad de corriente en x")
  # plt.xlabel("x")
  # plt.ylabel("y")
  # plt.imshow(jx[:,:],extent=[Xmin,Xmax,Ymin,Ymax],aspect='auto');plt.colorbar()

  # plt.subplot(1,2,2)
  # plt.title("Densidad de corriente en y")
  # plt.xlabel("x")
  # plt.ylabel("y")
  # plt.imshow(jy[:,:],extent=[Xmin,Xmax,Ymin,Ymax],aspect='auto');plt.colorbar()

  # # fig.savefig('my_figure.png')




"""# Resultados"""

plt.figure(1)
graficar_particula('iones', 50,'bo')
plt.figure(2)
graficar_particula('electrones', 50,'go')

plt.style.use('classic')


file_name = 'fields'+ str(i) +'.pkl'
df = pd.read_pickle(file_name)
print(file_name)
ex_tmp =  np.reshape(np.array(df['ex']),[Nx+1,Ny+1],order ='F')
ey_tmp =  np.reshape(np.array(df['ey']),[Nx+1,Ny+1], order ='F')
ez_tmp =  np.reshape(np.array(df['ez']),[Nx+1,Ny+1],order ='F')
bx_tmp =  np.reshape(np.array(df['bx']),[Nx+1,Ny+1],order ='F')
by_tmp =  np.reshape(np.array(df['by']),[Nx+1,Ny+1],order ='F')
bz_tmp =  np.reshape(np.array(df['bz']),[Nx+1,Ny+1],order ='F')
jx_tmp =  np.reshape(np.array(df['jx']),[Nx+1,Ny+1],order ='F')
jy_tmp =  np.reshape(np.array(df['jy']),[Nx+1,Ny+1],order ='F')
jz_tmp =  np.reshape(np.array(df['jz']),[Nx+1,Ny+1],order ='F')
divB_tmp = np.reshape(np.array(df['DivB']),[Nx+1,Ny+1],order ='F')


plt.figure(3);plt.clf()
plt.subplot(1,3,1)
plt.imshow(ex_tmp.T,extent=[Xmin,Xmax,Ymin,Ymax],aspect='auto', cmap='seismic');plt.colorbar()
plt.subplot(1,3,2)
plt.imshow(ey_tmp.T,extent=[Xmin,Xmax,Ymin,Ymax],aspect='auto', cmap='seismic');plt.colorbar()
plt.subplot(1,3,3)
plt.imshow(ez_tmp.T,extent=[Xmin,Xmax,Ymin,Ymax],aspect='auto', cmap='seismic');plt.colorbar()

plt.figure(4);plt.clf()
plt.subplot(1,3,1)
plt.imshow(bx_tmp.T,extent=[Xmin,Xmax,Ymin,Ymax],aspect='auto', cmap='seismic');plt.colorbar()
plt.subplot(1,3,2)
plt.imshow(by_tmp.T,extent=[Xmin,Xmax,Ymin,Ymax],aspect='auto', cmap='seismic');plt.colorbar()
plt.subplot(1,3,3)
plt.imshow(bz_tmp.T,extent=[Xmin,Xmax,Ymin,Ymax],aspect='auto', cmap='seismic');plt.colorbar()

plt.figure(5);plt.clf()
plt.subplot(1,3,1)
plt.imshow(jx_tmp.T,extent=[Xmin,Xmax,Ymin,Ymax],aspect='auto', cmap='seismic');plt.colorbar()
plt.subplot(1,3,2)
plt.imshow(jy_tmp.T,extent=[Xmin,Xmax,Ymin,Ymax],aspect='auto', cmap='seismic');plt.colorbar()
plt.subplot(1,3,3)
plt.imshow(jz_tmp.T,extent=[Xmin,Xmax,Ymin,Ymax],aspect='auto', cmap='seismic');plt.colorbar()


