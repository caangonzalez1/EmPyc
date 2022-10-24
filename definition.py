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

"""## CreaciC3n de arreglos"""

# crear malla
xx = np.linspace(Xmin, Xmax, Nx+1)
yy = np.linspace(Ymin, Ymax, Ny+1)
malla_x, malla_y = np.meshgrid(xx, yy)

energia_cinetica = np.zeros(n_iteraciones)
epsilon_r = np.zeros(n_iteraciones)

# Arreglos de Posiciones de electrones
x_BC = np.zeros([Np, dim])  # BC: Boundary condition position
x_act = np.zeros([Np, dim]) # act: actual (esta se usa para las condiciones iniciales)
x_new = np.zeros([Np, dim]) # new: nueva
# Arreglos de Velocidades
v_act = np.zeros([Np, dim])
v_new = np.zeros([Np, dim])

# Arreglos de Posiciones de iones
x_BC_i = np.zeros([Np, dim])  # BC: Boundary condition position
x_act_i = np.zeros([Np, dim]) # act: actual (esta se usa para las condiciones iniciales)
x_new_i = np.zeros([Np, dim]) # new: nueva
# Arreglos de Velocidades
v_act_i = np.zeros([Np, dim])
v_new_i = np.zeros([Np, dim])

# campo electromagnetico (B & E on the nodes)
bx = np.zeros([Nx+1, Ny+1])
by = np.zeros([Nx+1, Ny+1])
bz = np.zeros([Nx+1, Ny+1])
ex = np.zeros([Nx+1, Ny+1])
ey = np.zeros([Nx+1, Ny+1])
ez = np.zeros([Nx+1, Ny+1])
jx = np.zeros([Nx+1, Ny+1])
jy = np.zeros([Nx+1, Ny+1])
jz = np.zeros([Nx+1, Ny+1])
divB = np.zeros([Nx+1, Ny+1])

# Trabajo juan
bx_copy = np.zeros([Nx+1, Ny+1])
by_copy = np.zeros([Nx+1, Ny+1])
bz_copy = np.zeros([Nx+1, Ny+1])
ex_copy = np.zeros([Nx+1, Ny+1])
ey_copy = np.zeros([Nx+1, Ny+1])
ez_copy = np.zeros([Nx+1, Ny+1])

bx_tmp = np.zeros([Nx+1, Ny+1])
by_tmp = np.zeros([Nx+1, Ny+1])
bz_tmp = np.zeros([Nx+1, Ny+1])
ex_tmp = np.zeros([Nx+1, Ny+1])
ey_tmp = np.zeros([Nx+1, Ny+1])
ez_tmp = np.zeros([Nx+1, Ny+1])

# trabajo juan
xgc_bx=np.zeros([Nx,Ny])
xgc_ex=np.zeros([Nx,Ny])
xgc_ey=np.zeros([Nx,Ny])
xgc_by=np.zeros([Nx,Ny])
xgc_ez=np.zeros([Nx,Ny])
xgc_bz=np.zeros([Nx,Ny])
