import numpy as np
import pylab as pl
import scipy as sp
import bisect
from scipy.interpolate import interp1d
from scipy.interpolate import spline
from matplotlib import pyplot as plt

import pyfits

pl.rcParams['figure.figsize'] = (10.0, 7.0)
pl.rcParams['font.size'] = 20
pl.rcParams['font.family'] = 'serif'
pl.rcParams['lines.linewidth'] = 3

#MODIFICAR UBICACIONES
pathforfigs ='~/scripts/figscr'
pathforaux='/home/aguirre/aux'
filename=pathforaux+'/CascadeSpectra/Spectra/AtProduction_gammas.dat'
path=pathforaux+"/sensitivities/"


#MODIFICAR UBICACIONES DE LAS TABLAS
tabla50 = np.genfromtxt("/home/aguirre/scripts/tablasigmavmasaacabada.txt")
tabla502 = np.genfromtxt("/home/aguirre/scripts/simu2/tablasigmavmasa50.txt")
tabla = np.genfromtxt("/home/aguirre/scripts/simu1/tablasigmavmasa10.txt")
#tabla = np.genfromtxt("/home/aguirre/scripts/simu1/tablamedias10.txt")
tablaproc = np.genfromtxt("/home/aguirre/scripts/simu1/tablaproceeding.txt")
tablaprocscl = np.genfromtxt("/home/aguirre/scripts/simu1/tablasclW.txt")
tabla2 = np.genfromtxt("/home/aguirre/scripts/simu1/tablasigmavmasa11.txt")
tablaprocseg = np.genfromtxt("/home/aguirre/scripts/simu1/tablaseg1bb.txt")



Qe = 1.602176462e-19
TeV = 1
GeV = 1e-3 * TeV
MeV = 1e-6 * TeV
erg = 0.624151 * TeV 
eV = 1e-9 * GeV



def getTabla(tabla, tabla2=0,  index=1, medias = 0):
    #if medias == 0: #quitar esto cuando tenga todos los puntos
    tabla = tabla[tabla[:,0].argsort()]
    x = tabla[:,0]
    y = tabla[:,index]
    if medias == 1:
        tabla2 = tabla2[tabla2[:,0].argsort()] #con todos los puntos
        x = tabla2[:,0]
        aux = y
        y = tabla2[:,index]
        aux2 = [0 for i in xrange(len(y))]
        for i in xrange(len(y)):
            #y[i] = (aux[i]+y[i])/2  media lineal
            if aux[i] != y[i]:
                print aux[i]
                print y[i]
                aux2[i] = (aux[i]+y[i])/2
                y[i] = (y[i] - aux[i])/(np.log(y[i]) - np.log(aux[i]))
                print y[i]
                print aux2[i]
                print 'jeje'
    return x,y

    

 
fig=pl.figure(figsize=(15,10))

ax=fig.add_subplot(111)
ax.set_yscale('log')
ax.set_xscale('log')
#ax.set_xlim(0, 1)
#ax.set_ylim(0,1e6)
#ax.set_xlim(1e-5, 1)
#ax.set_ylim(1e-2,1e3)
ax.set_xlabel('$E$ [GeV]')
ax.set_ylabel('$<\sigma v>$ [cm$^3$ s$^{-1}$]')

ratioJsegscl = 0.14746

(E,sigmav) = getTabla(tabla50)
ax.plot(E, sigmav, label="datos 50h", color='orange', linewidth=1)
(E,sigmav) = getTabla(tabla502)
ax.plot(E, sigmav, label="datos 50h", color='#B84DFF', linewidth=1)
(E,sigmav) = getTabla(tabla)
ax.plot(E, sigmav, label="datos1 500h", color='red', linewidth=1)
#plotVarianza(ax,E,sigmav)
(E,sigmav) = getTabla(tabla2)
ax.plot(E, sigmav, label="datos2 500h", color='#00CCFF', linewidth=1)
#plotVarianza(ax,E,sigmav,fcolor='#20EEEE')
(E,sigmav) = getTabla(tablaproc)
ax.plot(E*1000, sigmav, label="Segue1 b", color='blue', linewidth=0.7,linestyle='--')
#(E,sigmav) = getTabla(tablaprocseg)
#ax.plot(E, sigmav, label="Segue1 b", color='pink', linewidth=1)
(E,sigmav) = getTabla(tablaprocscl)
ax.plot(E, sigmav*ratioJsegscl, label="Sculptor $W^+W^-$, ratio=0.147 (South)", color='black', linewidth=2)
#ax.plot(E, sigmav*0.33, label="SCL W prueba", color='#CC00CC', linewidth=1, linestyle='--')
(E,sigmav) = getTabla(tabla, tabla2, medias = 1)
ax.plot(E, sigmav, label="media datos 500h", color='#33CC33', linewidth=0.7, linestyle='--')
plt.legend(loc=4)
plt.show()
fig.savefig('sigmavsmi.png')
