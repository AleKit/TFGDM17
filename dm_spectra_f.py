import numpy as np
import pylab as pl
import scipy as sp
import bisect
from scipy.interpolate import interp1d
from scipy.interpolate import spline
from matplotlib import pyplot as plt

import pyfits

pl.rcParams['figure.figsize'] = (10.0, 7.0)
pl.rcParams['font.size'] = 18
pl.rcParams['font.family'] = 'serif'
pl.rcParams['lines.linewidth'] = 3


pathforfigs ='/home/ale/TFGF/'
pathforaux='/home/ale/TFGF'
filename=pathforaux+'/CascadeSpectra/Spectra/AtProduction_gammas.dat'
path=pathforaux+"/sensitivities/"


#vts_file = np.genfromtxt(path+"Instrument/VERITAS_V6_std_50hr_5sigma_VERITAS2014_DiffSens.dat")
vts_file = np.genfromtxt(path+"Instrument/VERITAS_ICRC2015_envelope.dat")

magic_file = np.genfromtxt(path+"Instrument/MAGIC_DiffSensCU.dat")

hess_file_combined = np.genfromtxt(path+"Instrument/HESS_August2015_CT15_Combined_Std.dat")
hess_file_stereo = np.genfromtxt(path+"Instrument/HESS_August2015_CT15_Stereo_Std.dat")
hess_file_envelope = np.genfromtxt(path+"Instrument/HESS_ICRC2015_envelope.dat")


hawc_1yr_file = np.genfromtxt(path+"Instrument/HAWC300_1y_QuarterDecade_DiffSens.dat")
hawc_5yr_file = np.genfromtxt(path+"Instrument/HAWC300_5y_QuarterDecade_DiffSens.dat")

fermi_b0_file = np.genfromtxt(path+"Instrument/fermi_lat_pass8_l0_b0.dat")
fermi_b30_file = np.genfromtxt(path+"Instrument/fermi_lat_pass8_l0_b30.dat")
fermi_b90_file = np.genfromtxt(path+"Instrument/fermi_lat_pass8_l0_b90.dat")

hs_file = np.genfromtxt(path+"Instrument/hiscore.dat")
lhaaso_file = np.genfromtxt(path+"Instrument/lhaaso.dat")


cta_n_file = np.genfromtxt(path+"North/CTA-Performance-North-50h-DiffSens.txt",skip_header=9)
cta_s_file = np.genfromtxt(path+"South/CTA-Performance-South-50h-DiffSens.txt",skip_header=9)


Qe = 1.602176462e-19
TeV = 1
GeV = 1e-3 * TeV
MeV = 1e-6 * TeV
erg = 0.624151 * TeV 
eV = 1e-9 * GeV

def Smooth(E, F):
    logE = np.log10(E)
    logF = np.log10(F)
    
    logEnew = np.linspace(logE.min(), logE.max(), 300)
    logF_smooth = spline(logE, logF, logEnew)
    
    Enew = [10**x for x in logEnew]
    F_smooth = [10**x for x in logF_smooth]
    
    return (Enew, F_smooth)

def getMAGIC(magic_file, Escale = GeV, smooth=True):
    x = 0.5*(magic_file[:,0] + magic_file[:,1]) * Escale
    y = magic_file[:,2]
    y_m = [y0 * 3.39e-11 * x0**(-2.51 - 0.21*np.log10(x0)) * x0 * x0 * 1e12 * Qe / 1e-7 for (x0, y0) in zip(x,y)]
    if smooth:
        return Smooth(x, y_m)
    else:
        return (x,y_m)

def getVERITAS(veritas_file, index=1, smooth=True):
    x = veritas_file[:,0]
    y = veritas_file[:,index]
    y_m = [y0 * 3.269e-11 * x0**(-2.474 - 0.191*np.log10(x0)) * x0 * x0 * 1e12 * Qe / 1e-7 for (x0, y0) in zip(x,y)]
    if smooth:
        return Smooth(x, y_m)
    else:
        return (x,y_m)
    

def getHESS(hess_file, smooth=True):
    x = hess_file[:,0]
    y = hess_file[:,1]
    
    x_m = [10**x0 for x0 in x]
    y_m = [y0 * x0 * x0 * 1e12 * Qe / 1e-7 for (x0,y0) in zip(x_m,y)]
    if smooth:
        return Smooth(x_m, y_m)
    else:
        return (x_m,y_m)
    
    
def getHESSEnvelope(hess_file, smooth=True):
    x = hess_file[:,0]
    y = hess_file[:,1]
    
    y_m = [y0 * x0 * x0 / erg for (x0,y0) in zip(x,y)]
    if smooth:
        return Smooth(x, y_m)
    else:
        return (x,y_m)    
    
def getHAWCFermi(hfile, Escale = GeV, smooth=True):
    # MeV scale for Fermi, GeV for HAWC
    x = hfile[:,0]
    y = hfile[:,1]
  
    if smooth:
        return Smooth(x * Escale, y)
    else:
        return (x * Escale,y)


def getCTA(ctafile, smooth=True):
    x = (ctafile[:,0] + ctafile[:,1])/2.
    y = ctafile[:,2]
    if smooth:
        return Smooth(x, y)
    else:
        return (x,y)

# Some useful units
#GeV = 1
#TeV = 1e3 * GeV
#erg = 624.15 * GeV
#eV = 1e-9 * GeV

def InterpolateTauEBL(E, redshift):
    #filename = "/a/home/tehanu/santander/ebl/ebl_z%0.1f.dat" % redshift
    filename = path + "ebl_z%0.1f.dat" % redshift
    eblfile = np.genfromtxt(filename)

    z = eblfile[:,0]
    Egamma = eblfile[:,1] * TeV
    Tau = eblfile[:,3]

    EBLInterp = interp1d(np.log10(Egamma), np.log10(Tau), kind='linear')

    TauValues = []

    for i in range(len(E)):
        if E[i] < Egamma[0]:
            TauValues.append(np.log10(Tau[0]))
        elif E[i] > Egamma[-1]:
            TauValues.append(np.log10(Tau[-1]))
        else:
            TauValues.append(EBLInterp(np.log10(E[i])))

    return [10**tau for tau in TauValues]

def SpectrumFlux(A, E, gamma, redshift = 0, Enorm = 1 * GeV, b = 0):
    if redshift > 0:
        tau = InterpolateTauEBL(E, redshift)
    else:
        tau = [0 for x in E]

    opacity = [np.exp(-t) for t in tau]

    return [A * (E0/Enorm)**(-gamma + b * np.log(E0/Enorm)) * exptau for (E0, exptau) in zip(E, opacity)]

def CrabSpectrumBroad(E):
    C = -0.12
    logf0 = -10.248
    Eic = 48*GeV # GeV
    a = 2.5
    return E**(-2)* 10**(logf0+C*(np.abs(np.log10(E/Eic)))**a)

def SpectrumIntegrator(Ec, Ewidths, Flux):
    nbins = 500
    Ebins = np.logspace(np.log10(Emin), np.log10(Emax), nbins)
    Ecenters = (Ebins[:-1]+Ebins[1:])/2
    Flux = SpectrumFlux(A, Ecenters, gamma)
    Ewidths = (Ebins[1:]-Ebins[:-1])
    return (Ecenters, np.sum(Ewidths * Flux))

def SpectrumIntAboveEnergy(Ec, Ewidths, Flux):
    prod = np.array([f * ew for (f, ew) in zip(Flux, Ewidths)])
    return [np.sum(prod[bisect.bisect(Ec,En):]) / TeV for En in Ec]

def plotSensitivity(ax, filename, legend="", xscale = 1, yscale=1., color='black', ls='-', lw=3, multE=False, delim=',', inCU=False, CrabEscale = 1):
    (Ec, flux) = plotCrabSpectrum(ax, scale=1, plot=False)
    f = interp1d(Ec, flux)

    if legend == "":
        legend = filename

    eblfile = np.genfromtxt(filename, delimiter=delim)

    x = eblfile[:,0] * xscale
    y = eblfile[:,1] * yscale

    l = zip(x,y)
    l.sort()

    xx = [x for (x,y) in l]
    yy = [y for (x,y) in l]

    if inCU:
        yy = f(xx) * yy * 1e-2

    if multE:
        zz = [xi * yi / TeV for (xi, yi) in zip(xx, yy)]
        ax.loglog(xx, zz, label=legend, linewidth=lw, color=color, ls=ls)
    else:
        ax.loglog(xx,yy, label=legend, linewidth=lw, color=color, ls=ls)    


e=1*GeV
e**2*CrabSpectrumBroad(e)


def plotIntegralSpectrum(ax, legend="", color='black', redshift=0, gamma=2, A=1e-8, Enorm = 1 * GeV, scale=1e-2, b = 0, fill=True, lwf=0.8, lwe=1, plot=True):
    Emin = 0.1 * GeV
    Emax = 1e8 * GeV
    nbins = 1500

    Ebins = np.logspace(np.log10(Emin), np.log10(Emax), nbins)
    Ec = (Ebins[:-1]+Ebins[1:])/2
    Ewidths = (Ebins[1:]-Ebins[:-1])
    Flux = SpectrumFlux(A, Ec, gamma, redshift, Enorm, b)
    IntFlux = SpectrumIntAboveEnergy(Ec, Ewidths, Flux)

    if fill:
        lowedge = [1e-16 for x in IntFlux]
        if plot:
            ax.fill_between(Ec,scale*Ec*IntFlux, lowedge, label="z = " + str(redshift),lw=0, alpha=0.08, color='#009933')
            ax.loglog(Ec,scale*Ec*IntFlux,lw=lwf, color='#009933', ls='-',alpha=0.5)
        return (Ec, scale*Ec*IntFlux)

    else:
        if plot:
            ax.loglog(Ec,scale*Ec*IntFlux,lw=lwe, color=color, ls='--')
        return (Ec, scale*Ec*IntFlux)


def plotSpectrum(ax, legend="", color='black', redshift=0, gamma=2, A=1e-8, Enorm = 1 * GeV, scale=1e-2, b = 0, fill=True, lwf=0.8, lwe=1, plot=True, fcolor='#009933', alpha=0.03):
    Emin = 0.1 * GeV
    Emax = 1e8 * GeV
    nbins = 1500

    Ebins = np.logspace(np.log10(Emin), np.log10(Emax), nbins)
    Ec = (Ebins[:-1]+Ebins[1:])/2
    Ewidths = (Ebins[1:]-Ebins[:-1])
    Flux = SpectrumFlux(A, Ec, gamma, redshift, Enorm, b)

    if fill:
        lowedge = [1e-16 for x in Flux]
        if plot:
            ax.fill_between(Ec,scale*Ec*Ec*Flux, lowedge, label="z = " + str(redshift),lw=0, alpha=alpha, color=fcolor)
            ax.loglog(Ec,scale*Ec*Ec*Flux,lw=lwf, color=color, ls='-',alpha=0.5)
        return (Ec, scale*Ec*Ec*Flux)

    else:
        if plot:
            ax.loglog(Ec,scale*Ec*Ec*Flux,lw=lwe, color=color, ls='--')
        return (Ec, scale*Ec*Ec*Flux)

def plotCrabSpectrumBroad(ax, legend="", color='black', scale=1, fill=True, lwf=0.8, lwe=1, plot=True, fcolor='grey', alpha=0.03):
    Emin = 0.1 * GeV
    Emax = 1e8 * GeV
    nbins = 1500

    Ebins = np.logspace(np.log10(Emin), np.log10(Emax), nbins)
    Ec = (Ebins[:-1]+Ebins[1:])/2
    Flux = CrabSpectrumBroad(Ec)
    
    if fill:
        lowedge = [1e-16 for x in Flux]
        if plot:
            ax.fill_between(Ec,scale*Ec*Ec*Flux, lowedge, lw=0, alpha=alpha, color=fcolor)
            ax.loglog(Ec,scale*Ec*Ec*Flux,lw=lwf, color=color, ls='-',alpha=0.5)
        return (Ec, scale*Ec*Ec*Flux)

    else:
        if plot:
            ax.loglog(Ec,scale*Ec*Ec*Flux,lw=lwe, color=color, ls='--')
        return (Ec, scale*Ec*Ec*Flux)

Ns=1e3
fullsky = 4 * np.pi


def getDMspectrum(option='e',finalstate='b',mass=1000,Jfactor=1.7e19,boost=1):
    #Options:
    #  e: outputs (E, dN/dE)
    #  e2: outputs (E, E**2 dN/dE)
    #  x: outputs (x,dN/dx)
    # mass in GeV
    # Jfactor in GeV2cm-5
    sigmav=3*1e-26 # annihilation cross section in cm3s-1
    data = np.genfromtxt (filename, names = True ,dtype = None,comments='#')

    massvals = data["mDM"]
    index = np.where(np.abs( (massvals - mass) / mass) < 1.e-3)
    xvals = 10**(data["Log10x"][index])
    
    
    def branchingratios(m_branon): #<sigmav>_particle / <sigmav>_total
    #PhysRevD.68.103505
        m_top = 172.44
        m_W = 80.4
        m_Z = 91.2
        m_h = 125.1
        m_c = 1.275
        m_b = 4.18
        m_tau = 1.7768
        if m_branon > m_top:
            c_0_top = 3.0 / 16 * m_branon ** 2 * m_top ** 2 * (m_branon ** 2 - m_top ** 2) * (1 - m_top ** 2 / m_branon ** 2) ** (1.0 / 2) 
        else:
            c_0_top = 0
        if m_branon > m_Z:
            c_0_Z = 1.0 / 64 * m_branon ** 2 * (1 - m_Z ** 2 / m_branon ** 2) ** (1.0 / 2) * (4 * m_branon ** 4 - 4 * m_branon ** 2 * m_Z ** 2 + 3 * m_Z ** 4)
        else:
            c_0_Z = 0
        if m_branon > m_W:
             c_0_W = 2.0 / 64 * m_branon ** 2 * (1 - m_W ** 2 / m_branon ** 2) ** (1.0 / 2) * (4 * m_branon ** 4 - 4 * m_branon ** 2 * m_W ** 2 + 3 * m_W ** 4)
        else:
             c_0_W = 0
        if m_branon > m_h:
             c_0_h = 1.0 / 64 * m_branon ** 2 * (2 * m_branon ** 2 + m_h ** 2) ** 2 * (1 - m_h ** 2 / m_branon ** 2) ** (1.0 / 2)
        else:
             c_0_h = 0
        if m_branon > m_c:
             c_0_c = 3.0 / 16 * m_branon ** 2 * m_c ** 2 * (m_branon ** 2 - m_c ** 2) * (1 - m_c ** 2 / m_branon ** 2) ** (1.0 / 2) 
        else:
             c_0_c = 0
        if m_branon > m_b:
             c_0_b = 3.0 / 16 * m_branon ** 2 * m_b ** 2 * (m_branon ** 2 - m_b ** 2) * (1 - m_b ** 2 / m_branon ** 2) ** (1.0 / 2) 
        else:
             c_0_b = 0
        if m_branon > m_tau:
             c_0_tau = 1.0 / 16 * m_branon ** 2 * m_tau ** 2 * (m_branon ** 2 - m_tau ** 2) * (1 - m_tau ** 2 / m_branon ** 2) ** (1.0 / 2) 
        else:
             c_0_tau = 0
        c_0_T = c_0_top + c_0_Z + c_0_W + c_0_h + c_0_c + c_0_b + c_0_tau
        br_t = (c_0_top / c_0_T)
        br_Z = c_0_Z / c_0_T
        br_W = c_0_W / c_0_T
        br_h = c_0_h / c_0_T
        br_c = c_0_c / c_0_T
        br_b = c_0_b / c_0_T
        br_tau = c_0_tau / c_0_T
        #f.append((c_0_T/(3*10**(-26)*math.pi**2))**(1./8))
        return {'masas': m_branon, 't': br_t, 'Z': br_Z, 'W': br_W, 'h': br_h, 'c': br_c, 'b': br_b, 'Tau': br_tau}
    
    #tau name modified in AtProduction_Gammas.dat

    if finalstate == "new":
        di = branchingratios(mass)
        flux_c = data[di.keys()[1]][index]/(np.log(10)*xvals) 
        flux_tau = data[di.keys()[2]][index]/(np.log(10)*xvals) 
        flux_b = data[di.keys()[3]][index]/(np.log(10)*xvals) 
        flux_t = data[di.keys()[4]][index]/(np.log(10)*xvals) 
        flux_W = data[di.keys()[5]][index]/(np.log(10)*xvals) 
        flux_Z = data[di.keys()[7]][index]/(np.log(10)*xvals) 
        flux_h = data[di.keys()[6]][index]/(np.log(10)*xvals) 
    
        loadspec_h = interp1d(xvals,flux_h)
        loadspec_Z = interp1d(xvals,flux_Z)
        loadspec_t = interp1d(xvals,flux_t)
        loadspec_W = interp1d(xvals,flux_W)
        loadspec_b = interp1d(xvals,flux_b)
        loadspec_c = interp1d(xvals,flux_c)
        loadspec_tau = interp1d(xvals,flux_tau)
    else:
        flux = data[finalstate][index]/(np.log(10)*xvals) #data is given in dN/d(log10(X)) = x ln10 dN/dx
        #flux = data[finalstate][index] 
        loadspec = interp1d(xvals,flux)

    def dNdx(x):
        fluxval = loadspec(x)
        if (x>1 or fluxval<0):
            return 0
        else:
            return fluxval
        
    def dNdx_new(x,di):
        fluxval_h = loadspec_h(x)
        if (x>1 or fluxval_h<0):
            fluxval_h = 0
        
        fluxval_Z = loadspec_Z(x)
        if (x>1 or fluxval_Z<0):
            fluxval_Z = 0
        
        fluxval_t = loadspec_t(x)
        if (x>1 or fluxval_t<0):
            fluxval_t = 0
        
        fluxval_W = loadspec_W(x)
        if (x>1 or fluxval_W<0):
            fluxval_W = 0
        
        fluxval_b = loadspec_b(x)
        if (x>1 or fluxval_b<0):
            fluxval_b = 0
        
        fluxval_c = loadspec_c(x)
        if (x>1 or fluxval_c<0):
            fluxval_c = 0
        
        fluxval_tau = loadspec_tau(x)
        if (x>1 or fluxval_tau<0):
            fluxval_tau = 0
        return (di.values()[1]*fluxval_c + di.values()[2]*fluxval_tau + 
                di.values()[3]*fluxval_b + di.values()[4]*fluxval_t +
               di.values()[5]*fluxval_W + di.values()[7]*fluxval_Z +
               di.values()[6]*fluxval_h)

    vdNdx = []
    x2vdNdx = []
    dNde = []
    e2dNde = []
    evals = []
    xvals2 = [] #aportacion mia
    if  option is 'e': #and boost > 1:
        #if mass == 5000:
        sigmavboost = sigmav * boost #no era necesario
        file1 = open("tabla"+str(mass)+str(finalstate)+str(sigmavboost)+".txt","w")

    logxvalsnew = np.linspace(-8.9,0,10000)
    xvalsnew = 10**logxvalsnew

    for i in range(len(xvalsnew)):
        x=xvalsnew[i]
        xvals2.append(x) #aportacion mia
        #vdNdx.append(dNdx(x))
        #x2vdNdx.append(x**2*dNdx(x))
        #dNde.append(dNdx(x)*Jfactor*GeV**2*sigmav*boost/(8*np.pi*(mass*GeV)**3))
        #e2dNde.append((1/erg)*x**2*dNdx(x)*Jfactor*GeV**2*sigmav*boost/(8*np.pi*mass*GeV))
        if finalstate == 'new':
            aux = dNdx_new(x,di)
        else:
            aux = dNdx(x)
        vdNdx.append(aux)
        x2vdNdx.append(x**2*aux)
        dNdeaux = aux*Jfactor*GeV**2*sigmav*boost/(8*np.pi*(mass*GeV)**3)
        dNde.append(dNdeaux)
        e2dNde.append((1/erg)*x**2*aux*Jfactor*GeV**2*sigmav*boost/(8*np.pi*mass*GeV))
        
        
        evals.append(x*mass*GeV)
        if option is 'e': #and boost > 1:
            #if mass == 5000 and dNdeaux != 0:
            if dNdeaux != 0:
                file1.write(str(x*mass*10**3) + " " + str(dNdeaux/(10**6)) + "\n")
                #print i
                #print(option, boost, mass, x*mass*10**3, dNdeaux/(10**6))
        #print(x, vdNdx[i], evals[i], e2dNde[i])
      #  if x == 1:
      #      break
    if option is 'e':
        #if mass == 5000 and boost > 1:
        file1.write(str(x*mass*10**3+1) + " " + "1e-99" + "\n")
        file1.write(str(x*mass*10**3+5) + " " + "1e-99" + "\n")
        file1.write(str(x*mass*10**3+10) + " " + "1e-99" + "\n")
        file1.close()
        return (evals,dNde)
    if option is 'e2':
        return (evals,e2dNde)
    if option is 'x':
        return (xvals2,vdNdx)
    if option is 'x2':
        return (xvals2,x2vdNdx)
    else:
        print('Option '+str(option)+' not supported')


fig=pl.figure(figsize=(15,10))

ax=fig.add_subplot(221)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlim(1e-7, 1)
ax.set_ylim(1e-7,1e6)
#ax.set_xlim(1e-5, 1)
#ax.set_ylim(1e-2,1e3)
ax.set_xlabel('$x$')
ax.set_ylabel('$dN/dx$')

(Edm,Fdm) = getDMspectrum('x','new',50)
ax.plot(Edm, Fdm, label="m = 0.05 TeV", color='red', linewidth=1)

(Edm,Fdm) = getDMspectrum('x','new',100)
ax.plot(Edm, Fdm, label="m = 0.1 TeV", color='blue', linewidth=1)

(Edm,Fdm) = getDMspectrum('x','new',150)
ax.plot(Edm, Fdm, label="m = 0.15 TeV", color='green', linewidth=1)

(Edm,Fdm) = getDMspectrum('x','new',250)
ax.plot(Edm, Fdm, label="m = 0.25 TeV", color='pink', linewidth=1)

(Edm,Fdm) = getDMspectrum('x','new',500)
ax.plot(Edm, Fdm, label="m = 0.5 TeV", color='#00CCFF', linewidth=1)

(Edm,Fdm) = getDMspectrum('x','new',1000)
ax.plot(Edm, Fdm, label="m = 1 TeV", color='#FF66FF', linewidth=1)

(Edm,Fdm) = getDMspectrum('x','new',5000)
ax.plot(Edm, Fdm, label="m = 5 TeV", color='#CC0066', linewidth=1)

(Edm,Fdm) = getDMspectrum('x','new',10000)
ax.plot(Edm, Fdm, label="m = 10 TeV", color='orange', linewidth=1)

(Edm,Fdm) = getDMspectrum('x','new',50000)
ax.plot(Edm, Fdm, label="m = 50 TeV", color='purple', linewidth=1)
plt.legend(loc=3, prop={'size':12}) 


ax=fig.add_subplot(223)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlim(1e-7, 1)
ax.set_ylim(1e-7,1)
ax.set_xlabel('$x$')
ax.set_ylabel('$x^2 dN/dx$')


#(Edm,Fdm) = getDMspectrum('x2','b',10000)
#ax.plot(Edm, Fdm, label="DM", color='pink', linewidth=1)

(Edm,Fdm) = getDMspectrum('x2','new',50)
ax.plot(Edm, Fdm, label="m = 0.05 TeV", color='red', linewidth=1)

(Edm,Fdm) = getDMspectrum('x2','new',100)
ax.plot(Edm, Fdm, label="m = 0.1 TeV", color='blue', linewidth=1)

(Edm,Fdm) = getDMspectrum('x2','new',150)
ax.plot(Edm, Fdm, label="m = 0.15 TeV", color='green', linewidth=1)

(Edm,Fdm) = getDMspectrum('x2','new',250)
ax.plot(Edm, Fdm, label="m = 0.25 TeV", color='pink', linewidth=1)

(Edm,Fdm) = getDMspectrum('x2','new',500)
ax.plot(Edm, Fdm, label="m = 0.5 TeV", color='#00CCFF', linewidth=1)

(Edm,Fdm) = getDMspectrum('x2','new',1000)
ax.plot(Edm, Fdm, label="m = 1 TeV", color='#FF66FF', linewidth=1)

(Edm,Fdm) = getDMspectrum('x2','new',5000)
ax.plot(Edm, Fdm, label="m = 5 TeV", color='#CC0066', linewidth=1)

(Edm,Fdm) = getDMspectrum('x2','new',10000)
ax.plot(Edm, Fdm, label="m = 10 TeV", color='orange', linewidth=1)

(Edm,Fdm) = getDMspectrum('x2','new',50000)
ax.plot(Edm, Fdm, label="m = 50 TeV", color='purple', linewidth=1)
plt.legend(loc=2, prop={'size':12}) 

ax=fig.add_subplot(222)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlim(2e-4, 60)
ax.set_ylim(5e-22,1e-5)
ax.set_xlabel('$E$ [TeV]')
ax.set_ylabel('$dN/dE$ [cm$^{-2}$ s$^{-1}$ TeV$^{-1}$]')

#(Edm,Fdm) = getDMspectrum('e','b',10)
#ax.plot(Edm, Fdm, label="DM", color='red', linewidth=1)

(Edm,Fdm) = getDMspectrum('e','new',50)
ax.plot(Edm, Fdm, label="m = 0.05 TeV", color='red', linewidth=1)

(Edm,Fdm) = getDMspectrum('e','new',100)
ax.plot(Edm, Fdm, label="m = 0.1 TeV", color='blue', linewidth=1)

(Edm,Fdm) = getDMspectrum('e','new',150)
ax.plot(Edm, Fdm, label="m = 0.15 TeV", color='green', linewidth=1)

(Edm,Fdm) = getDMspectrum('e','new',250)
ax.plot(Edm, Fdm, label="m = 0.25 TeV", color='pink', linewidth=1)

(Edm,Fdm) = getDMspectrum('e','new',500)
ax.plot(Edm, Fdm, label="m = 0.5 TeV", color='#00CCFF', linewidth=1)

(Edm,Fdm) = getDMspectrum('e','new',1000)
ax.plot(Edm, Fdm, label="m = 1 TeV", color='#FF66FF', linewidth=1)


(Edm,Fdm) = getDMspectrum('e','b',5000,boost=4e4) #######
(Edm,Fdm) = getDMspectrum('e','Tau',5000,boost=2e4) #######
(Edm,Fdm) = getDMspectrum('e','W',5000,boost=4e4) #######
(Edm,Fdm) = getDMspectrum('e','new',5000)
ax.plot(Edm, Fdm, label="m = 5 TeV", color='#CC0066', linewidth=1)


(Edm,Fdm) = getDMspectrum('e','new',10000)
ax.plot(Edm, Fdm, label="m = 10 TeV", color='orange', linewidth=1)

(Edm,Fdm) = getDMspectrum('e','new',50000)
ax.plot(Edm, Fdm, label="m = 50 TeV", color='purple', linewidth=1)
plt.legend(loc=3, prop={'size':12})

ax=fig.add_subplot(224)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlim(2e-4, 60)
ax.set_ylim(5e-22,1e-12)
ax.set_xlabel('$E$ [TeV]')
ax.set_ylabel('$E^2 dN/dE$ [erg cm$^{-2}$ s$^{-1}$]')

#(Edm,Fdm) = getDMspectrum('e2','b',10)
#ax.plot(Edm, Fdm, label="DM", color='red', linewidth=1)

(Edm,Fdm) = getDMspectrum('e2','new',50)
ax.plot(Edm, Fdm, label="m = 0.05 TeV", color='red', linewidth=1)

(Edm,Fdm) = getDMspectrum('e2','new',100)
ax.plot(Edm, Fdm, label="m = 0.1 TeV", color='blue', linewidth=1)

(Edm,Fdm) = getDMspectrum('e2','new',150)
ax.plot(Edm, Fdm, label="m = 0.15 TeV", color='green', linewidth=1)

(Edm,Fdm) = getDMspectrum('e2','new',250)
ax.plot(Edm, Fdm, label="m = 0.25 TeV", color='pink', linewidth=1)

(Edm,Fdm) = getDMspectrum('e2','new',500)
ax.plot(Edm, Fdm, label="m = 0.5 TeV", color='#00CCFF', linewidth=1)

(Edm,Fdm) = getDMspectrum('e2','new',1000)
ax.plot(Edm, Fdm, label="m = 1 TeV", color='#FF66FF', linewidth=1)

(Edm,Fdm) = getDMspectrum('e2','new',5000)
ax.plot(Edm, Fdm, label="m = 5 TeV", color='#CC0066', linewidth=1)

(Edm,Fdm) = getDMspectrum('e2','new',10000)
ax.plot(Edm, Fdm, label="m = 10 TeV", color='orange', linewidth=1)

(Edm,Fdm) = getDMspectrum('e2','new',50000)
ax.plot(Edm, Fdm, label="m = 50 TeV", color='purple', linewidth=1)

plt.legend(loc=3, prop={'size':12})
#plt.show()





fig=pl.figure()
ax=fig.add_subplot(111)
legends = []
ctain=True

(Emagic, Fmagic) = getMAGIC(magic_file)
ax.plot(Emagic, Fmagic, label="MAGIC", color='#A20025', linewidth=0.7)

#(Evts, Fvts) = getVERITAS(vts_file,index=1)
#vtsplot, = ax.plot(Evts, Fvts, label="VERITAS (50 hr)", color='red', linewidth=2)

#(Ehess, Fhess) = getHESS(hess_file_combined)
#ax.plot(Ehess, Fhess, label="HESS", color='#0050EF')

#(Ehess, Fhess) = getHESS(hess_file_stereo)
#ax.plot(Ehess, Fhess, label="HESS", color="#1BA1E2")

#(Ehess, Fhess) = getHESSEnvelope(hess_file_envelope)
#ax.plot(Ehess, Fhess, label="H.E.S.S.", color="#F0A30A",linewidth=4)

#(Ehawc, Fhawc) = getHAWCFermi(hawc_1yr_file)
#hawcplot1, = ax.plot(Ehawc, Fhawc, label="HAWC-300 - 1yr", color='#008A00', linewidth=2)

#(Ehawc, Fhawc) = getHAWCFermi(hawc_5yr_file)
#hawcplot2, = ax.plot(Ehawc, Fhawc, label="HAWC-300 - 5yr", color='#A4C400', linewidth=2)

#(Efermi, Ffermi) = getHAWCFermi(fermi_b0_file, Escale=MeV)
#ax.plot(Efermi, Ffermi, label="Fermi - b0", color='#bbbbbb')

#(Efermi, Ffermi) = getHAWCFermi(fermi_b30_file, Escale=MeV)
#ax.plot(Efermi, Ffermi, label="Fermi-LAT ($b = 30^{\circ})$ ", color='#00ABA9')

(Efermi, Ffermi) = getHAWCFermi(fermi_b90_file, Escale=MeV)
fermiplot, = ax.plot(Efermi, Ffermi, label="LAT (Pass8) - 10yr", color="#1BA1E2", linewidth=0.7)

#(Ehs, Fhs) = getHAWCFermi(hs_file, Escale=TeV)
#ax.plot(Ehs, Fhs, label="HiSCORE", color="#AA00FF", linestyle='-',linewidth=0.7)

#(Ehs, Fhs) = getHAWCFermi(lhaaso_file, Escale=TeV)
#ax.plot(Ehs, Fhs, label="LHAASO", color="#0050EF", linestyle='-',linewidth=0.7)

(Ecta, Fcta) = getCTA(cta_n_file)
ax.plot(Ecta, Fcta, label="CTA (50h, North)", linestyle='-', color='goldenrod',linewidth=1)

if ctain:
    (Ecta, Fcta) = getCTA(cta_s_file)
    ctaplot, = ax.plot(Ecta, Fcta, label="CTA (50h, South)", linestyle='-', color='#825A2C',linewidth=1)

#### Fermi IGRB ####

#figrb = np.genfromtxt(pathforaux+"/igrb.dat")

#Emean = 1000*(figrb[:,0] + figrb[:,1])/2.
#Ewidth = (figrb[:,1]-figrb[:,0]) * 1000

#Figrb = [4 * np.pi * (F/Ew) * scale * 1.60218e-6 * 1e-3 *  E**2 for (F, scale, E, Ew) in zip(figrb[:,2], figrb[:,5], Emean, Ewidth)]
#Figrb_err = [4 * np.pi * (F_err/Ew) * 1.60218e-6 * 1e-3 * scale * E**2 for (F_err, scale, E, Ew) in zip(figrb[:,3], figrb[:,5], Emean, Ewidth)]
#Figrb_err[-1] = 3e-14
#Flims = figrb[:,3] < 1e-3

#DNC ax.errorbar(Emean/1e6, Figrb, yerr=Figrb_err, xerr=Ewidth/2e6, marker='o',linestyle='',ecolor='red',color='red',mec='red',ms=3,uplims=Flims,capsize=3, linewidth=1)
#DNC ax.fill_between(Emean/1e6, [mean - err for (mean,err) in zip(Figrb, Figrb_err)], [mean + err for (mean,err) in zip(Figrb, Figrb_err)], zorder=0, alpha=0.5, color="#cccccc")
#ax.set_yscale('log')
#ax.set_xscale('log')
#ax.grid('on')

#if ctain:
#    first_legend = plt.legend(handles=[fermiplot,vtsplot,hawcplot1,hawcplot2,ctaplot], loc='upper left', bbox_to_anchor=(1.01, 1),fontsize=12)
#else:
#    first_legend = plt.legend(handles=[fermiplot,vtsplot,hawcplot1,hawcplot2], loc='upper left', bbox_to_anchor=(1.01, 1),fontsize=12)

#ax.add_artist(first_legend)
#legends.append(first_legend)
#ax.legend(bbox_to_anchor=(1.05, 1), ncol=1, loc=2, fontsize=14)

#ax.plot(hess_file[:,0], hess_file[:,1])
#ax.set_yscale('log')
#ax.set_xscale('log')

#ax.set_xlim(2e-4, 40)
#ax.set_ylim(4e-14,1e-10)
#ax.set_ylim(5e-14,1e-9)

#ax.set_xlabel('$E$ [TeV]')
#ax.set_ylabel('$E^2 d\Phi/dE$ [erg cm$^{-2}$ s$^{-1}$]')

#plotCrabSpectrum(ax, scale=1./(erg*GeV), fill=True, lwf=0.3, fcolor='#009933', alpha=0.05)
#plotCrabSpectrum(ax, scale=1e-1/(erg*GeV), fill=True, lwf=0.3, fcolor='#009933', alpha=0.05)
#plotCrabSpectrum(ax, scale=1e-2/(erg*GeV), fill=True, lwf=0.3, fcolor='#009933', alpha=0.05)
#plotCrabSpectrum(ax, scale=1e-3/(erg*GeV), fill=True, lwf=0.3, fcolor='#009933', alpha=0.05)


# Global fit nu
gamma=2.5
A0= (2/3.) * 6.7e-18 * fullsky / (erg*GeV)
Enorm = 100 * TeV

gamma=2.3
A0=1.5e-18 * fullsky / (GeV*erg)
Enorm = 100 * TeV

#gamma_wb = 2
#A_wb = 1e-8 * fullsky / GeV*erg
#Enorm_wb = 1 * GeV


#gamma=2.0
#A0=1.5e-18 * fullsky / (erg*GeV)
#Enorm = 100 * TeV

#plotSpectrum(ax, color='#009933',redshift=0, scale=1/Ns, gamma=gamma, A=A0, Enorm=Enorm, fcolor='#009933', alpha=0.1)
#plotSpectrum(ax, color='#009933',redshift=0.5, scale=1/Ns, gamma=gamma_wb, A=A_wb, Enorm=Enorm_wb, alpha=0.1, fcolor='#009933')
#plotSpectrum(ax, color='#666666',redshift=0.5, scale=1/Ns, gamma=gamma, A=A0, Enorm=Enorm)
#plotSpectrum(ax, color='#999999',redshift=1, scale=1/Ns, gamma=gamma, A=A0, Enorm=Enorm)



#plotSpectrum(ax, color='#000000',redshift=0, scale=1/Ns, gamma=gamma, A=A0, Enorm=Enorm)
#plotSpectrum(ax, color='#333333',redshift=0.1, scale=1/Ns, gamma=gamma, A=A0, Enorm=Enorm)
#plotSpectrum(ax, color='#666666',redshift=0.5, scale=1/Ns, gamma=gamma, A=A0, Enorm=Enorm)
#plotSpectrum(ax, color='#999999',redshift=1, scale=1/Ns, gamma=gamma, A=A0, Enorm=Enorm)
myalpha=0.015
myfcolor='blue'
mycolor='grey'
annotE=0.5*GeV  
myfontsize=10
myrot=30
if 1:
    plotCrabSpectrumBroad(ax,color=mycolor,scale=1/erg,alpha=myalpha,fcolor=myfcolor)
    ax.annotate('Crab', xy=(annotE,2*annotE**2*CrabSpectrumBroad(annotE)), xycoords='data',
                horizontalalignment='center', verticalalignment='center',fontsize=myfontsize,rotation=myrot)
if 1:
    plotCrabSpectrumBroad(ax,color=mycolor,scale=1e-1/erg,alpha=myalpha,fcolor=myfcolor)
    ax.annotate('10% Crab', xy=(annotE,2e-1*annotE**2*CrabSpectrumBroad(annotE)), xycoords='data',
                horizontalalignment='center', verticalalignment='center',fontsize=myfontsize,rotation=myrot)
if 1:
    plotCrabSpectrumBroad(ax,color=mycolor,scale=1e-2/erg,alpha=myalpha,fcolor=myfcolor)
    ax.annotate('1% Crab', xy=(annotE,2e-2*annotE**2*CrabSpectrumBroad(annotE)), xycoords='data',
                horizontalalignment='center', verticalalignment='center',fontsize=myfontsize,rotation=myrot)
if 1:
    plotCrabSpectrumBroad(ax,color=mycolor,scale=1e-3/erg,alpha=myalpha,fcolor=myfcolor)
    ax.annotate('0.1% Crab', xy=(annotE,2e-3*annotE**2*CrabSpectrumBroad(annotE)), xycoords='data',
                horizontalalignment='center', verticalalignment='center',fontsize=myfontsize,rotation=myrot)


if 1:

    hdulist1 = pyfits.open('spectrumdmbboost.fits')
    hdulist1.info()
    datos = hdulist1[1].data

    energ = datos['Energy']
    ed_energ = datos['ed_Energy']
    eu_energ = datos['eu_Energy']
    flux = datos['Flux']
    e_flux = datos['e_Flux']

    plt.errorbar(energ, flux, xerr=(ed_energ, eu_energ), yerr = e_flux, color = 'red', marker = 'o', label = r'spectrum $b\bar b$', fmt = '', zorder = 0)

    hdulist1.close()

    hdulist2 = pyfits.open('spectrumdmWboost.fits')
    hdulist2.info()
    datos = hdulist2[1].data

    energ = datos['Energy']
    ed_energ = datos['ed_Energy']
    eu_energ = datos['eu_Energy']
    flux = datos['Flux']
    e_flux = datos['e_Flux']

    plt.errorbar(energ, flux, xerr=(ed_energ, eu_energ), yerr = e_flux, color = 'green', marker = 'o', label = 'spectrum $W^+ W^-$', fmt = '', zorder = 0)

    hdulist2.close()

    mylinestyle='--'
    #mymass=3
    myboost=1
    (Edm,Fdm) = getDMspectrum('e2','b',5e3,boost=myboost)
    dmplot1 = ax.plot(Edm, Fdm, label=r"$m_\chi$ = "+str(5)+r" TeV ($b\bar b$, B$_f$=1e"+str("{:.1f}".format(np.log10(myboost)))+")", color='red', linewidth=1,linestyle=mylinestyle)

    (Edm,Fdm) = getDMspectrum('e2','Tau',5e3,boost=myboost)
    dmplot2 = ax.plot(Edm, Fdm, label=r"$m_\chi$ = "+str(5)+r" TeV ($\tau^- \tau^+$, B$_f$=1e"+str("{:.1f}".format(np.log10(myboost)))+")", color='blue', linewidth=1,linestyle=mylinestyle)

    (Edm,Fdm) = getDMspectrum('e2','W',5e3,boost=myboost)
    dmplot3 = ax.plot(Edm, Fdm, label=r"$m_\chi$ = "+str(5)+r" TeV ($W^+ W^-$, B$_f$=1e"+str("{:.1f}".format(np.log10(myboost)))+")", color='green', linewidth=1,linestyle=mylinestyle)
    myboost2 = 2e4
    myboost3 = 4e4
    (Edm,Fdm) = getDMspectrum('e2','b',5e3,boost=myboost3)
    dmplot4 = ax.plot(Edm, Fdm, label=r"$m_\chi$ = "+str(5)+r" TeV ($b\bar b$, B$_f$=1e"+str("{:.1f}".format(np.log10(myboost2)))+")", color='pink', linewidth=1,linestyle=mylinestyle)
    (Edm,Fdm) = getDMspectrum('e2','Tau',5e3,boost=myboost2)
    dmplot5 = ax.plot(Edm, Fdm, label=r"$m_\chi$ = "+str(5)+r" TeV ($\tau^- \tau^+$, B$_f$=1e"+str("{:.1f}".format(np.log10(myboost2)))+")", color='orange', linewidth=1,linestyle=mylinestyle)
    (Edm,Fdm) = getDMspectrum('e2','W',5e3,boost=myboost3)
    dmplot6 = ax.plot(Edm, Fdm, label=r"$m_\chi$ = "+str(5)+r" TeV ($W^+ W^-$, B$_f$=1e"+str("{:.1f}".format(np.log10(myboost3)))+")", color='purple', linewidth=1,linestyle=mylinestyle)

    plt.legend(loc=4, prop={'size':9}) #aportacion mia
    dmplots= []
    dmplots.append(dmplot1)
    dmplots.append(dmplot2)
    dmplots.append(dmplot3)
    dmplots.append(dmplot4)
    dmplots.append(dmplot5)
    dmplots.append(dmplot6)
    ax.set_xlim(1e-4, 1e3)
    ax.set_ylim(1e-16,1e-10)
    ax.set_xlabel('$E$ [TeV]')
    ax.set_ylabel('$E^2 dN/dE$ [erg cm$^{-2}$ s$^{-1}$]')
    #second_legend = plt.legend(handles=dmplots, ncol=1, loc='upper left',bbox_to_anchor=(1.01, .25), fontsize=12)
    #ax.add_artist(second_legend)
    #legends.append(second_legend)
#dummy = []
#aux_legend = plt.legend(handles=dummy,bbox_to_anchor=(1.5, .2),frameon=False)
#legends.append(aux_legend)
#fig.savefig(pathforfigs+'ic_sensitivities_cta_dm_'+str(mymass)+'_'+str(myboost)+'.pdf',bbox_extra_artists=legends,bbox_inches='tight')
plt.show() #aportacion mia
fig.savefig(pathforfigs+'ic_sensitivities_prueba.pdf',bbox_extra_artists=legends,bbox_inches='tight')
