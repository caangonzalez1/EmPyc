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

"""
## Velocidades y posiciones iniciales
"""
iz, d = 0, 0

# Ciclo para llenar posiciones y velocidades iniciales de los electrones y iones
for i in range(Np):
  # Posiciones en distribuciC3n uniforme de los iones y electrones
  x_act[i, 0]   = random.uniform(Xmin+iz, Xmax-d)
  x_act[i, 1]   = random.uniform(Xmin+iz, Xmax-d)
  x_act[i, 2]   = random.uniform(Xmin+iz, Xmax-d)
  x_act_i[i, 0] = random.uniform(Xmin+iz, Xmax-d)
  x_act_i[i, 1] = random.uniform(Xmin+iz, Xmax-d)
  x_act_i[i, 2] = random.uniform(Xmin+iz, Xmax-d)


for i in range(Np):
  v_act_i[i, 0] = random.gauss(0, 2*vthi)
  v_act_i[i, 1] = random.gauss(0, 2*vthi)
  v_act_i[i, 2] = random.gauss(0, 2*vthi)
  if(i<int(Np/2)):
    v_act[i, 0]   = random.gauss(-2*vthi, 2*vthe)
    v_act[i, 1]   = random.gauss(-2*vthi, 2*vthe)
    v_act[i, 2]   = random.gauss(-2*vthi, 2*vthe)
  else:
    v_act[i, 0]   = random.gauss(2*vthi, 2*vthe)
    v_act[i, 1]   = random.gauss(2*vthi, 2*vthe)
    v_act[i, 2]   = random.gauss(2*vthi, 2*vthe)

# Visualizar la distribuciC3n
#velDfi = pd.DataFrame(v_act_i)
#velDfe = pd.DataFrame(v_act)
#velDfi.plot.density()
#velDfe.plot.density()
# Visualizar la distribuciC3n
#velDfi = pd.DataFrame(x_act_i)
#velDfe = pd.DataFrame(x_act)
#velDfi.plot.density()
#velDfe.plot.density()


# Definicion de los campos uniformes
for j in range(0,Ny+1):
  for i in range(0,Nx+1):
    bz[i,j] = B_0
    
# Ciclo principal con campos uniformes!!!
t0 = 20
spread = 6






