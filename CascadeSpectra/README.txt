Supplementary material G. Elor, N. Rodd, T. Slatyer, W. Xue, 
“Model-Independent Indirect Detection Constraints on Hidden 
Sector Dark Matter” (2015)
Contact: Nick Rodd <nrodd@mit.edu>

The folder “spectra” contains files that can be used to load cascade
or direct spectra, the latter using the results of M. Cirelli et al., 
JCAP 1103, 051 (2011), 1012.4515.

In addition to this readme we provide two files with examples of how
to load the spectra:

  - CascadeSpectraExample.nb: how to load the spectra in Mathematica
  - CascadeSpectraExample.py: how to load the spectra in Python

Two notes on these files:
  1. .nb file was written in Mathematica 10. It runs fine (but with 
     some harmless error messages) in Mathematica 9. Due to the use 
     of the PlotLegends function, it will not operate with 
     Mathematica 8, but if these are deleted the files will run.
  2. .py file does not allow for interpolation to values of epsilon_f,
     m_phi or m_chi not in the various lists. If such a value is 
     required we suggest using the .nb file.

There are two filetypes in “spectra”:

===============================
======= Cascade Spectra =======
===============================
Number of files: 18
File naming: Cascade_<finalstate>_<spectratype>.dat
Description:
    Photon, positron and antiproton spectra for 7
    final states: gammas, e, mu, tau, b-quark, W-boson,
    Higgs and gluon (for antiproton only the last 4 of 
    these are provided).
Structure:
    Column [1]: value of Epsilon_f (or m_phi in
        the case of gluons or the positron spectrum 
	from gammas) - use care when going outside 
	the values of these parameters provided 
	which is 0.01, 0.03, 0.05, 0.07, 0.1, 0.2, 
	0.3, 0.4 and 0.5 for Epsilon_f and 10, 20, 
	40, 50, 80, 100, 500, 1000 and 2000 GeV for 
	m_phi.
    Column [2]: value of Log_10[x], where x=E/m_DM,
        running from -8.9 to 0 in 0.05 steps
    Column [3]: value of dN/dLog_10[x] for a 1-step
        cascade for that value of Epsilon_f (or 
	m_phi) and Log_10[x].
    Column [4-8]: as for column 3, but for 2- to 6-
        step cascade.

==============================
======= Direct Spectra =======
==============================
Number of files: 3
File naming: AtProduction_<spectratype>.dat
Description:
    Files created by M. Cirelli et al. that can be 
    used to extract the direct spectra. See that
    reference for details, or our example files for
    how to load.