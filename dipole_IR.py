#!/usr/bin/env python3
import sys
import numpy as np
from   scipy.io import netcdf
from   scipy.signal import correlate 

def load_charges(topfile, scale=1.):

    """
    Read charge of each atom in the system, units e.
    """
    
    f = open(topfile, "r")
    read = 0
    chargeList = []
    for line in f:
        if read == 1:
            read = 2
            continue
        if "FLAG CHARGE" in line:
            read = 1
            continue
        if read >= 2:
            if "%" in line:
                read = 0
            else:
                L = line.split()
                for s in L:
                    chargeList.append(float(s))
                
        if "FLAG ATOMIC_NUMBER" in line:
            read  = -1
            zList = []
            continue
        if read == -1:
            read = -2
            continue
        if read <= -2:
            if "%" in line:
                read = 0
                break
            L = line.split()
            for s in L:
                if s == "1":
                    zList.append(scale)
                else:
                    zList.append(1.0)
                
    return np.array(zList)*np.array(chargeList)/18.2223

def dipole_timeseries(chargeList, coordsNc):
    
    """
    Collect dipole moments time series from a netcdf file
    """
    f = netcdf.netcdf_file(coordsNc, 'r')
    n_atoms = len(chargeList)
    dipoles = []
    for frame in f.variables['coordinates']:
        coordSet  = np.array(frame)
        com       = np.sum(coordSet, axis=0)/n_atoms
        coordSet -= com
        px        = np.sum(coordSet[:,0] * chargeList)
        py        = np.sum(coordSet[:,1] * chargeList)
        pz        = np.sum(coordSet[:,2] * chargeList)
        dipoles.append(np.array([px,py,pz]))
    del frame
    return dipoles

def dipole_ft( dipoles  ):
    dipoles = np.array(dipoles)
    mod_dip = np.sqrt(np.sum(dipoles*dipoles, axis=1))
    print(np.shape(mod_dip))

    ft      = np.fft.fft(mod_dip)
    return ft

def acf(x):
    result = correlate(x, x, mode='same')
#    return result[len(result)//2:]
    return result
    
if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("require topfile, charge scale, timestep, and netcdf coordinates")
        raise ValueError
    top     = sys.argv[1]
    scale   = float(sys.argv[2])
    tstep   = float(sys.argv[3])

    charges = load_charges(top, scale)
    
    print("loaded %i charges" % len(charges))
    dipoles = []
    for i_nc in range(4,len(sys.argv)):
        dipoles += dipole_timeseries(charges, sys.argv[i_nc])
    dipoles = np.array(dipoles)
    for i in range(3):
        dipoles[:,i] -= np.mean(dipoles[:,i])

    
    n = len(dipoles)
    dip_acf      = np.zeros((n,4))
    dip_acf[:n//2,0] = np.arange(0,n//2)
    dip_acf[n//2:,0] = n - np.arange(n//2, n)
    
    dip_acf[:,1] = acf(dipoles[:,0])
    dip_acf[:,2] = acf(dipoles[:,1])
    dip_acf[:,3] = acf(dipoles[:,2])

    acf_fft      =  np.zeros((len(dip_acf), 5))

    
    ##0.1 fs timestep or smaller
    acf_fft[:,0] =  np.fft.fftfreq(len(dip_acf), tstep)
    acf_fft[:,1] =  acf_fft[:,0] / 3e10 ##wavenumber in cm^{-1}
    for i in range(1,4):
        ft =  np.fft.fft(dip_acf[:,i])
        acf_fft[:,i+1] = np.sqrt(ft.real*ft.real + ft.imag*ft.imag)
    N = len(acf_fft)
    np.savetxt("%s_IRcscale_%.2f.dat" % (sys.argv[-1],scale), acf_fft[:,:])
