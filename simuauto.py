import gammalib
import ctools
import os
import dm_spectra_f #borrar lo que no sean funciones
import random

#esta funcion calcula la seccion eficaz para la que la significancia de deteccion esta en torno a 25

def analisis(sigmav, mbranon,tobs,numtabla):
    file1 = open("tablasigmavmasa"+str(numtabla)+".txt","a")
    paso = 10.
    auxpaso = 2.3
    it = 1
    TS = 0 #inicializo
    while abs(TS-25) > 1:
#creamos la tabla con el espectro y su fichero XML
        dm_spectra_f.getDMspectrum("e","new",mbranon,boost=sigmav/(3*1e-26))
        os.system('mv tabla'+str(mbranon)+'new'+str(sigmav)+'.txt ~/scripts/tablas4may/tabla'+str(mbranon)+'sim'+str(it)+'n'+str(numtabla)+'.txt')
        nuevoarchivo = "$CTOOLS/share/models/dm"+str(mbranon)+'sim'+str(it)+'n'+str(numtabla)+".xml"
        straux = 's/ESCRIBIRAQUI/tabla' + str(mbranon) +'sim'+ str(it) +'n'+str(numtabla)+ "/"
        #print straux
        strsed = "sed '" + straux + "'" +' $CTOOLS/share/models/dmreemplazable.xml > ' + nuevoarchivo
        os.system(strsed)
#simulacion de ctools
        sim = ctools.ctobssim()
        sim["inmodel"] = nuevoarchivo 
        sim["outevents"] = "events.fits" #cambiar quiza
        sim["caldb"] = "prod2"
        sim["irf"] = "North_50h" #cambiar
        sim["ra"] = 10.11778
        sim["dec"] = 16.08194
        sim["rad"] = 5.0
        sim["tmin"] = 0.0
        sim["tmax"] = tobs # (param de entrada)
        sim["emin"] = 0.03
        sim["emax"] = 100.0
        sim["seed"] = int(100000*random.random())
        sim.execute()  #quiza mejor .run para no guardar estos archivos
        like = ctools.ctlike(sim.obs())
        like.run()
    # el valor de Test Statistic
        TSant = TS
        TS = like.obs().models()[0].ts()
        if TS > 25 and TSant < 25:
            paso = 1 + auxpaso
            auxpaso = auxpaso/2
            sigmav = sigmav/paso
        elif TS < 25 and TSant > 25:
            paso = 1 + auxpaso
            auxpaso = auxpaso/2
            sigmav = sigmav*paso
        elif TS < 25:
            sigmav = sigmav*paso
        elif TS > 25:
            sigmav = sigmav/paso
        it += 1
        if paso-1 < 1e-4:
            print ("paso muy peque")
            break
    file1.write(str(mbranon) + " " + str(sigmav) + "\n")
    file1.close()
    return sigmav
