import simuauto
import random
import numpy

#masa = input("masa:")
#sigmav = input("sigmav:")

#simuauto.analisis(sigmav,masa)
masas = numpy.logspace(1.7,5,20)

#numtabla = int(100000*random.random())
numtabla = 1 #cambiar a voluntad

for masa in masas:
    imasa = int(masa)
    simuauto.analisis(3e-24,imasa,1800000.0,numtabla) #cambiar a voluntad los parametros
