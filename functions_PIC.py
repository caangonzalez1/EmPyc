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

"""## PIC functions

### Velocity
"""
def velocity(x_act,v_act,ex,ey,ez,bx,by,bz,dt_tmp, q_s, m_s, number_particles):
    for j in range(number_particles):
        # Renombro posición y velocidades
        x0, y0, z0 = x_act[j, 0], x_act[j, 1], x_act[j, 2]
        ux0, uy0, uz0 = v_act[j, 0], v_act[j, 1], v_act[j, 2]
        # Interpolación: found the closest node indexes according to particle position
        ileft, iright, hxleft, jdown, jup, hydown = interpolation(x0, y0, Xmin, Ymin, dx, dy, xx, yy)
        # Obtener campos tiempo n
        Ex, Ey, Ez, Bx, By, Bz = get_fields_nodes(ex, ey, ez, bx, by, bz, ileft, iright, hxleft, jdown, jup, hydown)
        # BorisA
        v_new[j, 0],v_new[j, 1],v_new[j, 2] = BorisA(Ex, Ey, Ez, Bx, By, Bz, ux0, uy0, uz0, q_s, m_s, dt_tmp)
    return v_new

"""### Position"""

def position(x_act,v_new, x_new, x_BC, number_particles):
    for j in range(number_particles):
        x0, y0, z0 = x_act[j, 0], x_act[j, 1], x_act[j, 2]
        u_finalx, u_finaly, u_finalz = v_new[j, 0], v_new[j, 1], v_new[j, 2]
        x1, xBC1, x2, xBC2 , x3, xBC3 = leapfrog(x0, y0, z0, u_finalx, u_finaly, u_finalz)
        x_new[j, 0] = x1
        x_new[j, 1] = x2
        x_new[j, 2] = x3
        x_BC[j, 0] = xBC1
        x_BC[j, 1] = xBC2
        x_BC[j, 2] = xBC3
    
    return x_new,x_BC

"""### Energy"""

def energy(v_new):
    energia_cinetica = 0
    for j in range(Np):
        # Energía
        energia_cinetica = energia_cinetica + 0.5*m_e*(v_new[j,0]**2+v_new[j,1]**2+v_new[j,2]**2)
        #epsilon_r[i] = (energia_cinetica[i] - energia_cinetica[0])/energia_cinetica[0]
    energia_cinetica = energia_cinetica/Np
    return energia_cinetica

"""### Interpolation"""

def interpolation(posx, posy, xmin, ymin, dx, dy,xx, yy):
    ix = int((posx-xmin)/dx)  # index x
    xleft = posx  # position of the left node
    xright = xleft+dx  # position of the right node
    ileft = ix  # index of the left node
    iright = ix+1  # #index of the right node
    # particle fractional x-distance from the nearest node
    hxleft = (posx-xx[ix])/dx
    # y-direction
    jy = int((posy-ymin)/dy)  # index y
    ydown = posy  # position of the down node
    yup = ydown+dy  # position of the up node
    jdown = jy  # index of the down node
    jup = jy+1  # index of up node
    # particle fractional y-distance from the nearest node
    hydown = (posy-yy[jy])/dy
    # print(ileft,iright,hxleft,jdown,jup,hydown)
    return ileft, iright, hxleft, jdown, jup, hydown

"""### Get_fields_nodes"""

def get_fields_nodes(ex, ey, ez, bx, by, bz, ileft, iright, hxleft, jdown, jup, hydown):
    # weight functions
    w1 = (1.0-hxleft)*(1.0-hydown)
    w2 = hxleft*(1.0-hydown)
    w3 = (1.0-hxleft)*hydown
    w4 = hxleft*hydown
    # linear interpolations of fields on the particle position
    Ex = w1*ex[ileft, jdown] + w2*ex[iright, jdown] + w3*ex[ileft, jup] + w4*ex[iright, jup]
    Ey = w1*ey[ileft, jdown] + w2*ey[iright, jdown] + w3*ey[ileft, jup] + w4*ey[iright, jup]
    Ez = w1*ez[ileft, jdown] + w2*ez[iright, jdown] + w3*ez[ileft, jup] + w4*ez[iright, jup]
    Bx = w1*bx[ileft, jdown] + w2*bx[iright, jdown] + w3*bx[ileft, jup] + w4*bx[iright, jup]
    By = w1*by[ileft, jdown] + w2*by[iright, jdown] + w3*by[ileft, jup] + w4*by[iright, jup]
    Bz = w1*bz[ileft, jdown] + w2*bz[iright, jdown] + w3*bz[ileft, jup] + w4*bz[iright, jup]
    return Ex, Ey, Ez, Bx, By, Bz

"""### Boris"""

def BorisA(Ex, Ey, Ez, Bx, By, Bz, velx, vely, velz, q, m, dt):
    # Hallar la velocidad final
    # paso 1 ecuacion (3) del documento (Inicicializar velocidades)
    ux0 = velx
    uy0 = vely
    uz0 = velz
    # Velocidad un más adelante "(Ex*q/(2.0*m))""
    u_menosx = ux0 + (Ex*q/(2.0*m))
    u_menosy = uy0 + (Ey*q/(2.0*m))
    u_menosz = uz0 + (Ez*q/(2.0*m))
    # Ángulo de fase de la rotación
    # Paso 2 ecuacion (6) del documento
    tx = q*dt*Bx/(2.0*m)
    ty = q*dt*By/(2.0*m)
    tz = q*dt*Bz/(2.0*m)
    # paso 3 ecuacion (8) del documento.
    u_primax = u_menosx + ((u_menosy*tz)-(u_menosz*ty))
    u_primay = u_menosy + ((u_menosz*tx)-(u_menosx*tz))
    u_primaz = u_menosz + ((u_menosx*ty)-(u_menosy*tx))
    # Agrupando terminos:
    sx = 2*tx/(1.0+(tx*tx))
    sy = 2*ty/(1.0+(ty*ty))
    sz = 2*tz/(1.0+(tz*tz))
    # Paso 4 ecuacion(9)
    u_masx = u_menosx + ((u_primay*sz) - (u_primaz*sy))
    u_masy = u_menosy + ((u_primaz*sx) - (u_primax*sz))
    u_masz = u_menosz + ((u_primax*sy) - (u_primay*sx))
    # Paso 5
    u_finalx = u_masx + (Ex*q/(2*m))
    u_finaly = u_masy + (Ey*q/(2*m))
    u_finalz = u_masz + (Ez*q/(2*m))
    return u_finalx, u_finaly, u_finalz

"""### Leapfrog"""

def leapfrog(posx, posy, posz, ux, uy, uz):
    # Movimiento en x
    x_new = posx+ux*dt
    # Movimiento en y
    y_new = posy+uy*dt
    # Movimiento en z
    z_new = posz+uz*dt
    # Prueba para ver distancia excedida por dx o dy
    # if abs(x_new-posx)>dx:
    #   print(x_new-posx, 'cruzó más de una celda')
    # if abs(y_new-posy)>dy:
    #   print(y_new-posy, 'cruzó más de una celda')
    # Condiciones de frontera en x
    if x_new>Xmax:
        x_BC = x_new-Lx # x_BC: posición nueva con Boundary Condition
    elif x_new<Xmin:
        x_BC = x_new+Lx # x_BC: posición nueva con Boundary Condition
    else:
        x_BC = x_new

    # Condiciones de frontera en y
    if y_new>Ymax:
        y_BC = y_new-Ly # y_BC: posición nueva con Boundary Condition
    elif y_new<Ymin:
        y_BC = y_new+Ly # y_BC: posición nueva con Boundary Condition
    else:
        y_BC = y_new
        
    z_BC = z_new
    
    return x_new, x_BC, y_new, y_BC, z_new, z_BC

"""### Zig-Zag"""

def zig_zag(x_act, x_new, x_bc, ICx, ICy, peso_particula, q_s):
  for i in range(Np):
    x1, y1, z1 = x_act[i, 0], x_act[i, 1], x_act[i, 2]
    x2, y2, z2 = x_new[i, 0], x_new[i, 1], x_new[i, 2]
    x_BC, y_BC, z_BC = x_bc[i, 0], x_bc[i, 1], x_bc[i, 2]

    # print(x1, y1, x2, y2, x_BC, y_BC)

    ii0=int(math.floor((x1-Xmin)/dx))
    jj0=int(math.floor((y1-Ymin)/dy))
    ii1=int(math.floor((x2-Xmin)/dx))
    jj1=int(math.floor((y2-Ymin)/dy))
    ii2=int(math.floor((x_BC-Xmin)/dx))
    jj2=int(math.floor((y_BC-Ymin)/dy))

    # print(ii0,jj0,ii1,jj1,ii2,jj2)
    if ii2>64 or jj2>64:
      input()

    val1x=min(ii0*dx,ii1*dx)+dx;
    val2x=max(max(ii0*dx,ii1*dx),(x1+x2)/2.0);
    xr=min(val1x,val2x);

    val1y=min(jj0*dy,jj1*dy)+dy;
    val2y=max(max(jj0*dy,jj1*dy),(y1+y2)/2.0);
    yr=min(val1y,val2y);


#    if(ii0==ii1):
#      xr= 0.5*(x1+x2);
#    else:
#      xr= 0.5*((ii0+ii1)*dx);
        
#    if(jj0==jj1):
#      yr= 0.5*(y1+y2);
#    else:
#      yr= 0.5*((jj0+jj1)*dy);

    Fx1 = q_s*(xr-x1)/dt
    Fy1 = q_s*(yr-y1)/dt
    Fx2 = q_s*(x2-xr)/dt
    Fy2 = q_s*(y2-yr)/dt

    Wx1 = ((x1+xr)/2*dx) - ii0
    Wy1 = ((y1+yr)/2*dy) - jj0
    Wx2 = ((xr+x2)/2*dx) - ii1
    Wy2 = ((yr+y2)/2*dx) - jj1

    const = 1/(dx*dy)

    J1xU = const*Fx1*(1-Wy1)  # Jx(ii0+1/2, jj0)
    J1xD = const*Fx1*(Wy1)  # Jx(ii0+1/2, jj0+1)
    J1yL = const*Fy1*(1-Wx1)  # Jy(ii0, jj0+1/2)
    J1yR = const*Fy1*(Wx1)  # Jy(ii0+1, jj0+1/2)

    J2xU = const*Fx2*(1-Wy2)  # Jx(ii1+1/2, J2)
    J2xD = const*Fx2*(Wy2)  # Jx(ii1+1/2, J2+1)
    J2yL = const*Fy2*(1-Wx2)  # Jy(ii1, J2+1/2)
    J2yR = const*Fy2*(Wx2)  # Jy(ii1+1, J2+1/2)

    # Indices para interceldas que aporta la posición INICIAL

    ileft1=ii0
    iright1=ileft1+1
    jdown1=jj0
    jup1=jdown1+1

      # Llenar densidad de corriente en las interceldas (posición inicial)
    ICx[iright1][jdown1] = ICx[iright1][jdown1] + J1xU*peso_particula
    ICx[iright1][jup1]   = ICx[iright1][jup1]   + J1xD*peso_particula
    ICy[ileft1][jup1]    = ICy[ileft1][jup1]    + J1yL*peso_particula
    ICy[iright1][jup1]   = ICy[iright1][jup1]   + J1yR*peso_particula

    # Casos para los cuales la partícula se sale de la caja
    ileft2=ii2
    iright2=ileft2+1
    jdown2=jj2
    jup2=jdown2+1

    # Llenar densidad de corriente en las interceldas (posición final)
    
    # print(iright2, jdown2)
    ICx[iright2][jdown2] = ICx[iright2][jdown2] + J2xU*peso_particula
    ICx[iright2][jup2]   = ICx[iright2][jup2]   + J2xD*peso_particula
    ICy[ileft2][jup2]    = ICy[ileft2][jup2]    + J2yL*peso_particula
    ICy[iright2][jup2]   = ICy[iright2][jup2]   + J2yR*peso_particula



  for j in range(Ny):
    ICx[0,j] = ICx[Nx,j]
  for i in range(Nx):
    ICx[i,0]=ICx[i,Ny]+ICx[i,0]
    ICx[i,0]=ICx[i,Ny]

  for i in range(Nx):
    ICy[i,0] = ICy[i,Ny]
  for j in range(Ny):
    ICy[0,j]=ICy[Nx,j]+ICy[0,j]
    ICy[0,j]=ICy[Nx,j]

  # Normalización
  # ICx = ICx*I_0
  # ICy = ICy*I_0

  return ICx, ICy

