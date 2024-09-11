import numpy as np
import sncosmo
from scipy import special
import scipy
import random
from astropy.cosmology import Planck18
from astropy import units as u
import pandas as pd
import sncosmo
from astropy.table import Table
import matplotlib.pyplot as plt
import pickle
from os import listdir
from tqdm import trange
import ast

# Constants
C = 299792
SQ_ARCSECOND_ON_SKY = (41253) * (60**4)

def r_ein(sig, zl, zs):
    
    """
    Calculate the Einstein radius.

    Parameters:
    - sig (float): Lens mass dispersion.
    - zl (float): Lens redshift.
    - zs (float): Source redshift.

    Returns:
    - float: Einstein radius.
    """
    if zs<zl:
        return 0
    else: 
        ds=Planck18.angular_diameter_distance(zs).value
        dl=Planck18.angular_diameter_distance(zl).value
        dls=(ds*(1+zs)-dl*(1+zl))/(1+zs)
        rein=sig**2*((ds*C**2)/(206265*4*np.pi*dls))**-1
        return rein
        
def D_deltat(z_l,z_s):
"""
    Calculate time-delay distance.

    Parameters:
    - zl (float): Lens redshift.
    - zs (float): Source redshift.

    Returns:
    - float: Time-delay distance.
    """
    D_s=Planck18.angular_diameter_distance(z_s).value
    D_l=Planck18.angular_diameter_distance(z_l).value
    D_ls=(D_s*(1+z_s)-D_l*(1+z_l))/(1+z_s)
    Ddt=(1+z_l)*D_l*D_s/D_ls
    return Ddt
    
def lens_SN(lens, zs, catalogue_index, z_band):
    zl = lens[0]
    sigma = lens[1]

    Nsne = 1 #Number of positions per lens drawn

    R_e = r_ein(sigma, zl, zs)
    Ddt = D_deltat(zl, zs)
    td = Ddt * R_e**2 / 1000 #Catalogues are generated asusming Ddt of 1000 (kpc) and 

    q = round(lens[2], 2)

    if q > 0.99:
        q = 0.99
    if q < 0.2:
        q = 0.2

    df = files[q]
    df = df.sample(Nsne)
    df = df.reset_index()

    data_dict = {
        'zs':zs,
        'zl':zl,
        'RA_source' : (df['RA'].to_list()),
        'DEC_source' : (df['DEC'].to_list()),
        'tmax': (td * df['tmax']).tolist(),
        'max_sep': (R_e * df['max_sep']).tolist(),
        'mu_tot': df['mu_tot'].tolist(),
        'mu1': (np.array(df['mu_1'])).tolist(),
        't1': (td * np.array(df['t1'] - df['t1'])).tolist(),
        'mu2': df['mu_2'].tolist(),
        't2': (td * np.array(df['t2'] - df['t1'])).tolist(),
        'mu3': df['mu_3'].tolist(),
        't3': (td * np.array(df['t3'] - df['t1'])).tolist(),
        'mu4': df['mu_4'].tolist(),
        't4': (td * np.array(df['t4'] - df['t1'])).tolist(),
        'SN_index': catalogue_index,
        'peakMag_z': z_band
    }
    return data_dict
    
def process_catalog(catalog_path, z_band):
    catalogue = pd.read_csv(catalog_path)
    lenses = catalogue[['zl', 'sigma', 'ql']].values
    list_dicts = []  # List to hold dictionaries
    
    for idx in range(len(lenses)):
        lens_data = lens_SN(lenses[idx], catalogue['z_host'][idx], idx, z_band[idx])
        list_dicts.append(lens_data)
    
    return list_dicts  
    
files={}
qlist=[0.15]
for i in range(84):
    qlist.append(round(qlist[i]+0.01,2))

for q in qlist:
    files[q]=pd.read_csv(f'../Simulation_catalogues/positions_mu_t_q={q}.csv',index_col=0)
    
#########

catalog_path = '../outputs/SNIa/SNIa_GalaxyGalaxyPop.csv'
df=pd.read_csv(catalog_path)
mag_keys=['mag_g_transient']
for mag in mag_keys:
        df[mag] = [ast.literal_eval(val) for val in df[mag]]
df['peakMag_z'] = [np.nanmin(x) for x in df[mag_keys[0]]] 

DF_GLSN = process_bigger_catalog(catalog_path, df['peakMag_z']) 
DF_GLSN.to_csv('../outputs/SNIa/SNIa_larger_distribution.pkl'')
