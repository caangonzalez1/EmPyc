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

"""### Output: Archivos"""

def output_particle(x_act, v_act, i, especie):
  file_name = especie + str(i) + '.pkl'
  # Guardar posición y velocidades en pickle
  df_tmp = pd.DataFrame({'x':x_act[:, 0], 'y':x_act[:, 1], 'z':x_act[:, 2],
                          'vx':v_act[:, 0],'vy':v_act[:, 1],'vz':v_act[:, 2]})
  df_tmp.to_pickle(file_name) # convierte a .pickel
    
def output_fields(exp,eyp,ezp,bxp,byp,bzp,jxp,jyp,jzp,i,divBp, header = 'fields'):
    
  file_name = header+str(i)+'.pkl'
  ex_t = np.reshape(exp,[1, (Nx+1)*(Ny+1)])
  ey_t = np.reshape(eyp,[1, (Nx+1)*(Ny+1)])
  ez_t = np.reshape(ezp,[1, (Nx+1)*(Ny+1)])
  bx_t = np.reshape(bxp,[1, (Nx+1)*(Ny+1)])
  by_t = np.reshape(byp,[1, (Nx+1)*(Ny+1)])
  bz_t = np.reshape(bzp,[1, (Nx+1)*(Ny+1)])
  jx_t = np.reshape(jxp,[1, (Nx+1)*(Ny+1)])
  jy_t = np.reshape(jyp,[1, (Nx+1)*(Ny+1)])
  jz_t = np.reshape(jzp,[1, (Nx+1)*(Ny+1)])
  divB_t=np.reshape(divBp,[1, (Nx+1)*(Ny+1)])
  
  # Guardar posición y velocidades en pickle
  df_tmp = pd.DataFrame({'ex': ex_t[0,:], 'ey':ey_t[0,:], 'ez':ez_t[0,:],
                          'bx': bx_t[0,:], 'by':by_t[0,:], 'bz':bz_t[0,:],
                          'jx': jx_t[0,:], 'jy':jy_t[0,:], 'jz':jz_t[0,:],
                          'DivB':divB_t[0,:]})
  df_tmp.to_pickle(file_name)
  return

