from pynbody.analysis import profile
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
#import seaborn as sns
import numpy as np
import pynbody
import h5py

#from ..array import SimArray
from pynbody import config
import logging
logger = logging.getLogger('pynbody.analysis.ionfrac')

#from scipy.interpolate import interp3d
from scipy.interpolate import RegularGridInterpolator as rgi

from pynbody.analysis.interpolate import _interpolate3d

def interpolate3d(x, y, z, x_vals, y_vals, z_vals, vals):
    """
    Interpolate on a 3D regular grid. 
    Yields results identical to scipy.interpolate.interpn. 
    Input
    -----
    x,y,z : points where the interpolation will be performed
    x_vals, y_vals, z_vals : xyz values of the reference grid
    vals : grid values
    """

    # cast x_vals, y_vals and z_vals to float64

    x_vals = x_vals.astype(np.float64)
    y_vals = y_vals.astype(np.float64)
    z_vals = z_vals.astype(np.float64)
    vals = vals.astype(np.float64)

    result_array = np.empty(len(x), dtype=np.float64)

    _interpolate3d.interpolate3d(len(x),
                                 x, y, z,
                                 len(x_vals), x_vals,
                                 len(y_vals), y_vals,
                                 len(z_vals), z_vals,
                                 vals,
                                 result_array)

    return result_array

def N_OVI(f):
    ovi = pynbody.analysis.ionfrac_edit.calculate(f,ion='ovi',mode='new')
    m_p = 1.6726 * 10**-24 #g
    #print(ovi)
    #print(f.gas['OxMassFrac'])
    #print(f.gas['rho'].in_units('g cm**-3')*ovi*f.gas['OxMassFrac']/(16*m_p))
    return f.gas['rho'].in_units('g cm**-3')*ovi*f.gas['OxMassFrac']/(16*m_p)

from pynbody.analysis.interpolate import interpolate3d
import h5py

def hdf5_ion_frac(sim, element, ion):
    """hdf5_ion_frac(sim, element, ion)')
       'sim=hg, element="C", ion=3 (for CIV)')"""
    #if ion == 'oi':
    #    print('Loading OI')
    #    nion = 0
        
    iffile = './hm2012_hr.h5'
    ifs = h5py.File(iffile)
    
    x_vals = ifs[element].attrs['Parameter1']   # HYDROGEN Densities
    y_vals = ifs[element].attrs['Parameter2']   # Redshifts
    z_vals = ifs[element].attrs['Temperature']  # Temperatures

    vals = ifs[element][ion]
    y = np.zeros(len(sim.gas))
    x = np.log10(sim.gas['rho'].in_units('m_p cm^-3')).view(np.ndarray)
    y[:] = sim.properties['z']
    y = y
    z = np.log10(sim.gas['temp']).view(np.ndarray)
 
    n = len(sim.gas)
    n_z_vals = len(z_vals)
    n_y_vals = len(y_vals)
    n_x_vals = len(x_vals)
    
    # get values off grid to minmax                 
    x[np.where(x < np.min(x_vals))] = np.min(x_vals)
    x[np.where(x > np.max(x_vals))] = np.max(x_vals)
    y[np.where(y < np.min(y_vals))] = np.min(y_vals)
    y[np.where(y > np.max(y_vals))] = np.max(y_vals)
    z[np.where(z < np.min(z_vals))] = np.min(z_vals)
    z[np.where(z > np.max(z_vals))] = np.max(z_vals)
    
    # interpolate                              
    result_array = interpolate3d(x, y, z, x_vals, y_vals, z_vals, vals)

    return 10 ** result_array


def hdf5_ion_frac_justox(sim, ion):
    if ion == 'oi':
        print('Loading OI')
        nion = 0
    elif ion == 'oii':
        print('Loading OII')
        nion = 1
    elif ion == 'oiii':
        print('Loading OIII')
        nion = 2
    elif ion == 'oiv':
        print('Loading OIV')
        nion = 3
    elif ion == 'ov':
        print('Loading OV')
        nion = 4
    elif ion == 'ovi':
        print('Loading OVI')
        nion = 5
    elif ion == 'ovii':
        print('Loading OVII')
        nion = 6
    elif ion == 'oviii':
        print('Loading OVIII')
        nion = 7
    else:
        print('Specified ion incompatible; Try again.')
        
    iffile = '/nobackup/nnsanche/hm2012_hr.h5'
    ifs = h5py.File(iffile)
    
    x_vals = ifs['O'].attrs['Parameter1']   # HYDROGEN Densities
    y_vals = ifs['O'].attrs['Parameter2']   # Redshifts
    z_vals = ifs['O'].attrs['Temperature']  # Temperatures

    vals = ifs['O'][nion]
    y = np.zeros(len(sim.gas))
    x = np.log10(sim.gas['rho'].in_units('m_p cm^-3')).view(np.ndarray)
    y[:] = sim.properties['z']
    y = y
    z = np.log10(sim.gas['temp']).view(np.ndarray)
 
    n = len(sim.gas)
    n_z_vals = len(z_vals)
    n_y_vals = len(y_vals)
    n_x_vals = len(x_vals)
    
    # get values off grid to minmax                 
    x[np.where(x < np.min(x_vals))] = np.min(x_vals)
    x[np.where(x > np.max(x_vals))] = np.max(x_vals)
    y[np.where(y < np.min(y_vals))] = np.min(y_vals)
    y[np.where(y > np.max(y_vals))] = np.max(y_vals)
    z[np.where(z < np.min(z_vals))] = np.min(z_vals)
    z[np.where(z > np.max(z_vals))] = np.max(z_vals)
    
    # interpolate                              
    result_array = interpolate3d(x, y, z, x_vals, y_vals, z_vals, vals)


    return 10 ** result_array



if __name__ == '__main__':

    k = 0
    sim = ['/nobackupp2/nnsanche/pioneer50h243.1536g1bwK1BH/pioneer50h243.1536gst1bwK1BH.004096','/nobackupp2/nnsanche/pioneer50h243GM1.1536gst1bwK1BH/pioneer50h243GM1.1536gst1bwK1BH.004096','/nobackupp2/nnsanche/pioneer50h243GM4.1536gst1bwK1BH/pioneer50h243GM4.1536gst1bwK1BH.004096','/nobackup/nnsanche/pioneer50h243GM5.1536gst1bwK1BH/pioneer50h243GM5.1536gst1bwK1BH.004096','/nobackupp2/nnsanche/pioneer50h243GM6.1536gst1bwK1BH/pioneer50h243GM6.1536gst1bwK1BH.004096','/nobackup/nnsanche/pioneer50h243GM7.1536gst1bwK1BH/pioneer50h243GM7.1536gst1bwK1BH.004096']
    labels = ['P0','GM1','GM4','GM5','GM6','GM7']
    colors = sns.cubehelix_palette(8)
    Z_sun = 0.0142 # (Asplund 2009; https://arxiv.org/pdf/0909.0948.pdf)
    m_p = 1.6726 * 10**-24 #g
    print('LOADING SIM:',labels[k])
    
    ######################
    # READ IN SIMULATION #
    ######################
    f = pynbody.load(sim[k])
    pynbody.analysis.halo.center(f.star)
    f.physical_units()
    h = f.halos()
    h1 = h[1]
    pynbody.analysis.angmom.faceon(h1)
    
    ###################
    # ISOLATE CGM GAS #   
    ###################
    # Isolate and remove disk stars within radius 0-10 kpc & vertically 10 kpc 
    r_max = 10  # kpc
    twenty_kpc_incm = 6.171*(10**22)
    #z_max = 10 #4 # kpc
    
    Rg_d = ((h1.g['x'].in_units('kpc'))**2. + (h1.g['y'].in_units('kpc'))**2. + (h1.g['z'].in_units('kpc'))**2.)**(0.5)
    disk_gas_xyzmax =  (Rg_d < r_max)
    #disk_gas_zmax  = (h1.g['z'].in_units('kpc') < z_max) & (h1.g['z'].in_units('kpc') > -z_max)
    disk_gas_mask = disk_gas_xyzmax #& disk_gas_zmax
    disk_gas = h1.g[disk_gas_mask] #& disk_gas_zmax]
    CGM_gas  = h1.g[~disk_gas_mask]
    
    CGM_gas['ovi_npz'] = pynbody.analysis.ionfrac.calculate(CGM_gas,ion='ovi',mode='new') 
    CGM_gas['ovi_hdf5'] = hdf5_ion_frac(CGM_gas,ion='ovi')
#    quit()
    
    CGMprofile = profile.Profile(CGM_gas,min='0.1 kpc',max='250 kpc')
    
    # Column Densities
    Novi_npz = np.log10((CGMprofile['mass'].in_units('g')*CGMprofile['OxMassFrac']*CGMprofile['ovi_npz']/(16*m_p))/CGMprofile._binsize.in_units('cm**2'))
    Novi_hdf5 = np.log10((CGMprofile['mass'].in_units('g')*CGMprofile['OxMassFrac']*CGMprofile['ovi_hdf5']/(16*m_p))/CGMprofile._binsize.in_units('cm**2'))
    
    print('USING .NPZ DATA:')
    #print('Total mass in OVI in CGM:', np.sum(CGM_gas['OxMassFrac']*CGM_gas['mass']*CGM_gas['ovi_npz']))
    print('OVI fractions:',np.average(CGM_gas['ovi_npz']))
    
    print('USING HDF5 DATA:')
    #print('Total mass in OVI in CGM:', np.sum(CGM_gas['OxMassFrac']*CGM_gas['mass']*ovi))
    print('OVI fractions:',np.average(CGM_gas['ovi_hdf5']))
    
    plt.plot(CGMprofile['rbins'].in_units('kpc'),Novi_hdf5,color='Salmon')
    plt.plot(CGMprofile['rbins'].in_units('kpc'),Novi_npz,color='SteelBlue')
    
    plt.title('z = 0.17')
    plt.ylabel(r'log(N$_{ovi}$) [cm$^{-2}$]',size=15)
    plt.xlabel('R [kpc]',size=15)
    plt.ylim(12.5,17.5)
    plt.xlim(-10,260)
    plt.legend(ncol=2,fontsize=15)
    #plt.savefig('ALLGMs_plusnoBH_Novi_R.pdf')
    plt.show()
    plt.close()
    
    
    quit()
    OVI = CGM_gas['rho'].in_units('g cm**-3')*ovi*CGM_gas['OxMassFrac']/(16*m_p)
    print('Total mass in CGM:', np.sum(CGM_gas['mass']))
    print('Total mass in metals:',np.sum(CGM_gas['mass']*CGM_gas['metals'])) # 'metals' *IS* metallicity
    print('Total mass Oxygen in CGM:', np.sum(CGM_gas['OxMassFrac']*CGM_gas['mass']),CGM_gas['mass'].units)
    print('Total mass in OVI in CGM:', np.sum(CGM_gas['OxMassFrac']*CGM_gas['mass']*ovi))
    #print('OVI Density: ',OVI,OVI.units)
    print('OVI fractions:',np.average(ovi))
    
    quit()
