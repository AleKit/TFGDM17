#------------------------------------------------
# HEADER
#------------------------------------------------
##
## CASCADE SPECTRA
##
## This file describes the contents of the "Spectra" directory and further 
## below examples of how to load them into Python
##
## Supplementary material for G. Elor, N. Rodd, T. Slatyer, and W. Xue,
## "Model-Independent Indirect Detection Constraints on Hidden Sector Dark
## Matter" (2015)
## Contact: Nick Rodd <nrodd@mit.edu>
##
#------------------------------------------------
##
## Overview
##
## Files within "Spectra" contain the \[Gamma], e^+/e^- and p spectrum of 
## dark matter (DM) annihilations leading to 1-6 step cascades ending in 
## e, mu, tau, b-quark, W, h or gluon final states. A separate .dat file 
## is provided for each type of spectrum and final state. Note all spectra 
## here are calculated using the large hierarchy approximation, which can 
## break down in certain situations. The contents of these files are 
## designed to be similar to the results of M. Cirelli et al., JCAP 1103, 
## 051 (2011), 1012.4515.
##
#------------------------------------------------
##
## Content of files
##
## There are 18 files of the form Cascade_{Final State}_{Spectrum Type}.dat. 
## Each contains 8 columns and 1612 rows (the first row contains the column 
## labels and all others contain numerical values). For all files except 
## gluon final states and the positron spectrum from gammas, the columns are: 
## 1. epsilon_f value; 
## 2. Log10(x), where x=E/Subscript[m, DM]; 
## 3-8. the value of dN/dLog10(x)=ln(10)*x*dN/dx of an n=1 (column 3) up to 
## n=6 (column 8) spectrum at that value of epsilon_f and x. 
## For final state gluons and the positron spectrum from gammas, the columns 
## are:
## 1. m_phi value;
## 2. Log10(x), where x=E/Subscript[m, DM]; 
## 3-8. the value of dN/dLog10(x)=ln(10)*x*dN/dx of an n=1 (column 3) up to 
## n=6 (column 8) spectrum at that value of m_phi and x. 
## The reason for this difference in the case of gluons and the positron 
## spectrum from gammas is that we have to use m_phi values instead of 
## epsilon_f=2m_f/m_phi as for these final states m_f=0 making epsilon_f 
## less useful.
##
#------------------------------------------------
##
## Range of Parameters
##
## Log10(x) ranges from -8.9 to 0 in steps on 0.05, covering the entire range 
## where x^2 dN/dx for each spectra is non-zero. epsilon_f is evaluated at 0.01, 
## 0.03, 0.05, 0.07, 0.1, 0.2, 0.3, 0.4 and 0.5, whilst m_phi is given at 10, 20, 
## 40, 50, 80, 100, 500, 1000 and 2000 GeV.
##
## NB: if using a value of epsilon_f or m_phi not in this list - especially if it
## is outside the range of values provided - it is recommended that linear 
## interpolation be used.
##
## NB: as we are using an interpolating function to determine the spectrum, it can
## occasionally go negative. This is unphysical can cause issues if the Log of the 
## spectrum is being used, and so we set the spectrum to 0 wherever the 
## interpolating function would send it negative.
##
#------------------------------------------------
##
## Brief description of files
##
## Photon Spectra:
## Cascade_{Gam,E,Mu,Tau,B,W,H,G}_gammas.dat - photon spectrum from final state
## {photons, electrons, muons, taus, b-quarks, Ws, Higgs, gluons}
##
## Positron/Electron Spectra:
## Cascade_{Gam,E,Mu,Tau,B,W,H,G}_positrons.dat - positron (or equivalently 
## electron) spectrum from final state {photons, electrons, muons, taus, b-quarks, 
## Ws, Higgs, gluons}
##
## Antiproton Spectra:
## Cascade_{B,W,H,G}_antiprotons.dat - antiproton spectrum from final state {b-quarks,
## Ws, Higgs, gluons}
## NB: electron, muon and tau antiproton spectra are negligible and so not given
##
## Direct Spectra:
## AtProduction_{gammas,positrons,antiprotons}.dat - Direct {photon, positron, 
## antiproton} spectrum as calculated by M. Cirelli et al. These files were produced 
## by those authors, but we provide them here for convenience.
##
#------------------------------------------------
# ENDHEADER
#------------------------------------------------

# Import required modules
from scipy.interpolate import interp1d
import numpy as np
from matplotlib import pyplot as plt

#------------------------------------------------
# Example 1
# Photon Cascade Spectra of Final State W Bosons for epsilon_f=0.1
#------------------------------------------------
##
## In this example we load the photon spectrums of 1-6 step cascade DM 
## annihilations terminating in W bosons for epsilon_f=0.1.

filename = './Spectra/Cascade_W_gammas.dat'

with open(filename) as f:
    lines = (line for line in f if not line.startswith('#'))
    data = np.genfromtxt (lines, names = True ,dtype = None)

epsvals = data["EpsF"]
eps_f = 0.1
index = np.where(np.abs( (epsvals - eps_f) / eps_f) < 1.e-3)
xvals = 10**(data["Log10x"][index])

flux = [data["n"+str(i)][index]/(np.log(10)*xvals) for i in range(1,7)]
loadspec = [interp1d(xvals,flux[i]) for i in range(6)]
def dNdx(x,step):
    fluxval = loadspec[step-1](x)
    if (x>1 or fluxval<0):
        return 0
    else:
        return fluxval

plt.plot(xvals,[x**2*dNdx(x,1) for x in xvals],label='n=1',color='Purple')
plt.plot(xvals,[x**2*dNdx(x,2) for x in xvals],label='n=2',color='Blue')
plt.plot(xvals,[x**2*dNdx(x,3) for x in xvals],label='n=3',color='Green')
plt.plot(xvals,[x**2*dNdx(x,4) for x in xvals],label='n=4',color='Pink')
plt.plot(xvals,[x**2*dNdx(x,5) for x in xvals],label='n=5',color='Orange')
plt.plot(xvals,[x**2*dNdx(x,6) for x in xvals],label='n=6',color='Red')
plt.title('Ex. 1: Photon Cascade Spectra into W Bosons for $\epsilon_f=0.1$',fontsize=14)
plt.xscale('log')
plt.xlabel('$x$', fontsize=18)
plt.ylabel('$x^2 dN / dx$', fontsize=18)
plt.ylim([0.0,0.12])
plt.xlim([10**-6,1])
plt.legend(fontsize=12,loc=2)
plt.show()

#------------------------------------------------
# Example 2
# Antiproton Cascade Spectra of Final State Gluons for m_phi=50 GeV
#------------------------------------------------
##
## In this example we load the antiproton spectrums of 1-6 step cascade DM 
## annihilations terminating in Gluons for m_phi=50 GeV. 

filename = './Spectra/Cascade_G_antiprotons.dat'

with open(filename) as f:
    lines = (line for line in f if not line.startswith('#'))
    data = np.genfromtxt (lines, names = True ,dtype = None)

mphivals = data["mPhiGeV"]
mphi = 50
index = np.where(np.abs( (mphivals - mphi) / mphi) < 1.e-3)
xvals = 10**(data["Log10x"][index])

flux = [data["n"+str(i)][index]/(np.log(10)*xvals) for i in range(1,7)]
loadspec = [interp1d(xvals,flux[i]) for i in range(6)]
def dNdx(x,step):
    fluxval = loadspec[step-1](x)
    if (x>1 or fluxval<0):
        return 0
    else:
        return fluxval

plt.plot(xvals,[x**2*dNdx(x,1) for x in xvals],label='n=1',color='Purple')
plt.plot(xvals,[x**2*dNdx(x,2) for x in xvals],label='n=2',color='Blue')
plt.plot(xvals,[x**2*dNdx(x,3) for x in xvals],label='n=3',color='Green')
plt.plot(xvals,[x**2*dNdx(x,4) for x in xvals],label='n=4',color='Pink')
plt.plot(xvals,[x**2*dNdx(x,5) for x in xvals],label='n=5',color='Orange')
plt.plot(xvals,[x**2*dNdx(x,6) for x in xvals],label='n=6',color='Red')
plt.title('Ex. 2: Antiproton Cascade Spectra into Gluons for $m_{\phi}=50$ GeV',fontsize=14)
plt.xscale('log')
plt.xlabel('$x$', fontsize=18)
plt.ylabel('$x^2 dN / dx$', fontsize=18)
plt.ylim([0.0,0.04])
plt.xlim([10**-6,1])
plt.legend(fontsize=12,loc=2)
plt.show()


#------------------------------------------------
# Example 3
# Positron Direct Spectrum of Final State b-quarks for m_chi=100 GeV
#------------------------------------------------
##
## In this example we load the 0-step positron spectrum for a DM 
## annihilation proceeding directly  into b-quarks for epsilon_f=0.3.
## For other final states, replace the b in finalstate = "b" with
## e (electrons), Mu (muons), Tau (taus), b (b-quarks), W (W bosons),
## g (gluons), Gamma (photons) or h (Higgs 

filename = './Spectra/AtProduction_positrons.dat'

finalstate = "b"

with open(filename) as f:
    lines = (line for line in f if not line.startswith('#'))
    data = np.genfromtxt (lines, names = True ,dtype = None)

massvals = data["mDM"]
mass = 100
index = np.where(np.abs( (massvals - mass) / mass) < 1.e-3)
xvals = 10**(data["Log10x"][index])

flux = data[finalstate][index]/(np.log(10)*xvals)
loadspec = interp1d(xvals,flux)
def dNdx(x):
    fluxval = loadspec(x)
    if (x>1 or fluxval<0):
        return 0
    else:
        return fluxval

plt.plot(xvals,[x**2*dNdx(x) for x in xvals],label='n=0',color='dimgray')
plt.title('Ex. 3: Positron Direct Spectrum into b-quarks for $m_{\chi}=100$ GeV',fontsize=14)
plt.xscale('log')
plt.xlabel('$x$', fontsize=18)
plt.ylabel('$x^2 dN / dx$', fontsize=18)
plt.ylim([0.0,0.06])
plt.xlim([10**-6,1])
plt.legend(fontsize=12,loc=2)
plt.show()
