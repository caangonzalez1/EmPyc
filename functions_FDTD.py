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

from parameters import *
from definition import *

"""## FDTD functions

### FDTD_Halftime
"""
def FDTD_halftime(ex,ey,ez,bx,by,bz,dt_tmp):
    for i in range(Nx+1):
        for j in range(1,Ny+1):
            bx[i,j]=bx[i,j]-((dt_tmp/dy)*(ez[i,j]-ez[i,j-1]))
    #Condiciones de frontera Bx: Periodicidad
    for i in range(Nx+1):
        bx[i,0]=bx[i,Ny]

    for i in range(1,Nx+1):
        for i in range(Ny+1):
            by[i,j]=by[i,j]+((dt_tmp/dx)*((ez[i,j]-ez[i-1,j])))
    #Condiciones de frontera By: Periodicidad
    for j in range(Ny+1):
        by[0,j]=by[Nx,j]

    for i in range(Nx):
        for j in range(Ny):
            bz[i,j]= bz[i,j] + ( ((dt_tmp/dy)*((ex[i,j+1]-ex[i,j]))) - ((dt_tmp/dx)*(ey[i+1,j]-ey[i,j])) )
    #Condiciones de frontera Bz: Periodicidad
    #Esquina
    bz[Nx,0]=bz[0,0]
    for i in range(Nx+1):
        bz[i,Ny]=bz[i,0]
    for j in range(Ny+1):
        bz[Nx,j]=bz[0,j]

    return bx,by,bz

"""### FDTD"""
def FDTD(ex,ey,ez,bx,by,bz,jx,jy,jz,ex_copy,ey_copy,ez_copy,bx_copy,by_copy,bz_copy):

    #Evolución del campo BX,BY Y BZ ----> dt
    #Evolución de Bx,By,Bz ---> dt

    for i in range(Nx+1):
      for j in range(1,Ny+1):
        ex_copy[i,j],ey_copy[i,j],ez_copy[i,j]=ex[i,j],ey[i,j],ez[i,j]
        bx_copy[i,j],by_copy[i,j],bz_copy[i,j]=bx[i,j],by[i,j],bz[i,j]


    for i in range(Nx+1):
        for j in range(1,Ny+1):
            bx[i,j]=bx[i,j] - ((dt/dy)*(ez[i,j]-ez[i,j-1]))
    #Condiciones de frontera Bx: Periodicidad
    for i in range(Nx+1):
        bx[i,0]=bx[i,Ny]

    for i in range(1,Nx+1):
        for j in range(Ny+1):
            by[i,j]=by[i,j] + ((dt/dx)*((ez[i,j]-ez[i-1,j])))
    #Condiciones de frontera By: Periodicidad
    for j in range(Ny+1):
        by[0,j]=by[Nx,j]

    for i in range(Nx):
        for j in range(Ny):
            bz[i,j]= bz[i,j] + ( ((dt/dy)*((ex[i,j+1]-ex[i,j]))) - ((dt/dx)*(ey[i+1,j]-ey[i,j])) )
    #Condiciones de frontera Bz: Periodicidad
    #Esquina
    bz[Nx,0]=bz[0,0]
    for i in range(Nx+1):
        bz[i,Ny]=bz[i,0]
    for j in range(Ny+1):
        bz[Nx,j]=bz[0,j]

    #Evolución de Ex,Ey,Ez ------> dt

    for i in range(Nx+1):
        for j in range(1,Ny+1):
            ex[i,j]= ex[i,j] + ((dt/dy)*(bz[i,j]-bz[i,j-1])) - jx[i,j]
    #Condiciones de frontera Bx: Periodicidad
    for i in range(Nx+1):
        ex[i,0]=ex[i,Ny]

    for i in range(1,Nx+1):
        for j in range(Ny+1):
            ey[i,j]= ey[i,j] - ((dt/dx)*((bz[i,j]-bz[i-1,j]))) - jy[i,j]
    #Condiciones de frontera By: Periodicidad
    for j in range(Ny+1):
        ey[0,j]=ey[Nx,j]

    for i in range(Nx):
        for j in range(Ny):
            ez[i,j]=ez[i,j] - ( ((dt/dy)*((bx[i,j+1]-bx[i,j]))) - ((dt/dx)*(by[i+1,j]-by[i,j])) ) - jz[i,j]
    #Condiciones de frontera Bz: Periodicidad
    #Esquina
    ez[Nx,0]=ez[0,0]
    for i in range(Nx+1):
        ez[i,Ny]=ez[i,0]
    for j in range(Ny+1):
        ez[Nx,j]=ez[0,j]

    for i in range(Nx+1):
        for j in range(Ny+1):
          bx_copy[i,j]=0.5*(bx_copy[i,j]+bx[i,j])
          by_copy[i,j]=0.5*(by_copy[i,j]+by[i,j])
          bz_copy[i,j]=0.5*(bz_copy[i,j]+bz[i,j])

    return ex,ey,ez,bx,by,bz,ex_copy,ey_copy,ez_copy,bx_copy,by_copy,bz_copy

"""### Divergence"""

def divBt(bx_tmp,by_tmp,bz_tmp,divB,i):
  divB = np.zeros([Nx+1,Ny+1])
  for ii in range(1,Nx+1):
    for jj in range(1,Ny+1):
      #print(i,j)
      divB[ii,jj] = (bx_tmp[ii,jj] - bx_tmp[ii-1,jj]) + (by_tmp[ii,jj] - by_tmp[ii,jj-1])
  
  print('t=',i*dt,'<devB>=',np.sqrt(np.mean(divB**2)))
  
  return divB

"""### Interc&Nodes Functions"""

#lleva los valores de los campos en los nodes(nx_1,ny+1) a los centros  (nx,ny)
def nodes_to_c(xgc_ex,xgc_ey,xgc_ez,xgc_bx,xgc_by,xgc_bz,bx,by,ex,ey,bz,ez):
  for i in range(Nx):
    for j in range(Ny):
      xcg_bx=0.25*(bx[i,j]+bx[i+1,j]+bx[i,i+1]+bx[i+1,j+1])
      xcg_ex=0.25*(ex[i,j]+ex[i+1,j]+ex[i,i+1]+ex[i+1,j+1])
      xcg_ey=0.25*(ey[i,j]+ey[i+1,j]+ey[i,i+1]+ey[i+1,j+1])
      xcg_by=0.25*(by[i,j]+by[i+1,j]+by[i,i+1]+by[i+1,j+1])
      xcg_bx=0.25*(bx[i,j]+bx[i+1,j]+bx[i,i+1]+bx[i+1,j+1])
      xcg_ez=0.25*(ez[i,j]+ez[i+1,j]+ez[i,i+1]+ez[i+1,j+1])
      xcg_bz=bz[i,j]


#Lleva los valores de los campos en los centros (nx,ny) a los nodes nodos
def c_to_nodes(xcg_ex, xcg_ey, xcg_ez, xcg_bx,xcg_by,xcg_bz,bx,by,ex,ey,bz,ez):
  for i in range(1,Nx):
    for j in range(1,Ny):
      ex[i,j]=0.25*(xgc_ex[i-1,j-1]+xgc_ex[i,j-1]+xgc_ex[i-1,j]+xgc_ex[i,j])
      bx[i,j]=0.25*(xgc_bx[i-1,j-1]+xgc_bx[i,j-1]+xgc_bx[i-1,j]+xgc_bx[i,j])
      by[i,j]=0.25*(xgc_by[i-1,j-1]+xgc_by[i,j-1]+xgc_by[i-1,j]+xgc_by[i,j])
      ey[i,j]=0.25*(xgc_ey[i-1,j-1]+xgc_ey[i,j-1]+xgc_ey[i-1,j]+xgc_ey[i,j])
      ez[i,j]=0.25*(xgc_ez[i-1,j-1]+xgc_ez[i,j-1]+xgc_ez[i-1,j]+xgc_ez[i,j])
    
  #for Ez 2D setup
  for i in range(1,Nx):
    for j in range(1,Ny):
      bz[i,j]=xgc_bz[i-1,j-1]

  for i in range(Nx):
    ex[i,Ny]=ex[i,0]
    bx[i,Ny]=bx[i,0]
    ey[i,Ny]=ey[i,0]
    by[i,Ny]=by[i,0]
    ez[i,Ny]=ez[i,0]

  for j in range(Ny):
    ex[Nx,j]=ex[0,j]
    bx[Nx,j]=bx[0,j]
    ey[Nx,j]=ey[0,j]
    by[Nx,j]=by[0,j]
    ez[Nx,j]=ez[0,j]

  for i in range(Nx+1):
    bz[i,0]=bz[i,Ny]

  for j in range(Ny+1):
    bz[0,j]=bz[Nx,j]

  #Esquina
  bz[0,0]=bz[Nx,Ny]


  for i in range(Nx):
    for j in range(Ny):
      bx_copy[i,j]=bx[i,j]
      by_copy[i,j]=by[i,j]
      bz_copy[i,j]=bz[i,j]
      ex_copy[i,j]=ex[i,j]
      by_copy[i,j]=ey[i,j]
      ez_copy[i,j]=ez[i,j]

  
def nodes_to_inter(ICx,ICy,bx_copy,by_copy,bz_copy,ex_copy,ey_copy,ez_copy,ex,ey,ez,bx,by,bz):
  for i in range(Nx):
    for j in range(Ny):
      bx_copy[i,j]=bx[i,j]
      by_copy[i,j]=by[i,j]
      bz_copy[i,j]=bz[i,j]
      ex_copy[i,j]=ex[i,j]
      by_copy[i,j]=ey[i,j]
      ez_copy[i,j]=ez[i,j]

  for i in range(Nx+1):
    for j in range(Ny+1):
      ey_copy[i,j]=0.5*(ey[i,j]+ey[i-1,j])
      by_copy[i,j]=0.5*(by[i,j]+by[i-1,j])
      ICx[i,j]=0.5*(ICx[i,j]+ICx[i,j-1])

  #Condiciones de frontera

  for j in range(Ny+1):
    ey[0,j]=ey[Nx,0]
    by[0,j]=by[Nx,0]
    ICy[0,j]=ICy[Nx,0]

    for i in range(Nx+1):
      for j in range(Ny+1):
        ex_copy[i,j]=0.5*(ex[i,j]+ex[i-1,j])
        bx_copy[i,j]=0.5*(bx[i,j]+bx[i-1,j])
        ICy[i,j]=0.5*(ICy[i,j]+ICy[i,j-1])

    for i in range(Nx+1):
      ex[i,0]=ey[i,Ny]
      bx[i,0]=by[i,Ny]
      ICx[i,0]=ICx[i,Ny]

def interc_to_nodes(ex_copy,ey_copy,ez_copy, bx_copy,by_copy,bz_copy, ex_nodes, ey_nodes, ez_nodes, bx_nodes, by_nodes, bz_nodes):
  
  for i in range(Nx+1):
    for j in range(Ny):
      ex_nodes[i,j]=0.5*(ex_copy[i,j]+ex_copy[i,j+1])
      bx_nodes[i,j]=0.5*(bx_copy[i,j]+bx_copy[i,j+1])

  for i in range(Nx+1):
    ex_nodes[i,0]=ex_copy[i,Ny]
    bx_nodes[i,0]=bx_copy[i,Ny]

  for i in range(Nx):
    for j in range(Ny+1):
      ey_nodes[i,j]=0.5*(ey_copy[i,j]+ey_copy[i+1,j])
      by_nodes[i,j]=0.5*(by_copy[i,j]+by_copy[i+1,j])

  for j in range(Ny+1):
    ey_nodes[0,j]=ex_copy[Nx,j]
    by_nodes[0,j]=bx_copy[Nx,j]

  ez_nodes, bz_nodes =  ez_copy,bz_copy

  return ex_nodes, ey_nodes, ez_nodes, bx_nodes, by_nodes, bz_nodes


def energy_fields(bx_tmp,by_tmp,bz_tmp):
    energia_f = 0.5*(np.mean(bx_tmp**2 + bx_tmp**2 + bz_tmp**2 ))
    return energia_f
