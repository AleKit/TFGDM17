import pylab as pl
from matplotlib import pyplot as plt
import math


def branchingratios(m_branon_i, m_branon_f): #<sigmav>_particle / <sigmav>_total
    #PhysRevD.68.103505
    m_top = 172.44
    m_W = 80.4
    m_Z = 91.2
    m_h = 125.1
    m_c = 1.275
    m_b = 4.18
    m_tau = 1.7768
    br_t=[]
    br_Z=[]
    br_W=[]
    br_h=[]
    br_c=[]
    br_b=[]
    br_tau=[]
    masas = []
    f = []
    for m_branon in xrange(m_branon_i,m_branon_f):
        masas.append(m_branon)
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
        br_t.append(c_0_top / c_0_T)
        br_Z.append(c_0_Z / c_0_T)
        br_W.append(c_0_W / c_0_T)
        br_h.append(c_0_h / c_0_T)
        br_c.append(c_0_c / c_0_T)
        br_b.append(c_0_b / c_0_T)
        br_tau.append(c_0_tau / c_0_T)
        f.append((c_0_T/(3* 10**4 / (197.3**2 * 299792458)*math.pi**2))**(1./8))
    return {'masas': masas, 't': br_t, 'Z': br_Z, 'W': br_W, 'h': br_h, 'c': br_c, 'b': br_b, 'tau': br_tau, 'f': f}


fig=pl.figure(figsize=(14,6))
pl.rcParams['font.size'] = 17
ax=fig.add_subplot(121)

ax.set_yscale('log')
ax.set_xscale('log')
#ax.set_xlim(0, 1)
#ax.set_ylim(0,1e6)
#ax.set_xlim(1e-5, 1)
#ax.set_ylim(1e-2,1e3)
ax.set_xlabel('mass (GeV)')
ax.set_ylabel('Br')

brs = branchingratios(2,5000)
#print brs.keys()
ax.plot(brs.values()[0], brs.values()[1], label="c", color='red', linewidth=1)
ax.plot(brs.values()[0], brs.values()[2], label="$\\tau$", color='blue', linewidth=1)
ax.plot(brs.values()[0], brs.values()[3], label="b", color='green', linewidth=1)
ax.plot(brs.values()[0], brs.values()[4], label="t", color='pink', linewidth=1)
ax.plot(brs.values()[0], brs.values()[5], label="W", color='yellow', linewidth=1)
ax.plot(brs.values()[0], brs.values()[7], label="h", color='purple', linewidth=1)
ax.plot(brs.values()[0], brs.values()[8], label="Z", color='orange', linewidth=1)

plt.legend()

ax=fig.add_subplot(122)
ax.set_yscale('log')
ax.set_xscale('log')
#ax.set_xlim(0, 1)
#ax.set_ylim(0,1e6)
#ax.set_xlim(1e-5, 1)
#ax.set_ylim(1e-2,1e3)
ax.set_xlabel('mass (GeV)')
ax.set_ylabel('f (GeV)')
ax.plot(brs.values()[0], brs.values()[6], label="f", color='black', linewidth=1)

plt.show()
