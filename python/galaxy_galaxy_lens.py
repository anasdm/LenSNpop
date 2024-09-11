# Load necessary libraries
import numpy as np
import pandas as pd
from FastLensSim import FastLensSim
from tqdm import trange
import random
from astropy.cosmology import Planck18

# Load DataFrame, can be output file from combine_lens_and_source.py
df = pd.read_csv('../outputs/SNIa/SNIa_GalaxyGalaxyPop.csv') 

# Initialize arrays for computed results
magnification_gal = np.zeros(len(df))
x = np.zeros(len(df))
y = np.zeros(len(df))

# Initialize FastLensSim instance outside the loop
S = FastLensSim('ideal')
S.bfac = 2
S.rfac = 2

for i in trange(len(df)):
    bs = df['R_Ein'][i]
    ml = {'g_SDSS': df['mag_g_SDSS'][i], 'r_SDSS': df['mag_r_SDSS'][i], 'i_SDSS': df['mag_i_SDSS'][i]}
    rl = {'g_SDSS': df['rl'][i], 'r_SDSS': df['rl'][i], 'i_SDSS': df['rl'][i]}
    ql = df['ql'][i]
    ms = {'g_SDSS': df['mag_g_host'][i], 'r_SDSS': df['mag_r_host'][i], 'i_SDSS': df['mag_i_host'][i]}
    qs = (1 - df['e'][i]) / (1 + df['e'][i])
    ps = df['a_rot'][i]
    rs = df['rs'][i]

    xs = df['R_Ein'][i] * df['RA_SN'][i] + (random.choice([1, -1])) * df['ra_off'][i]
    ys = df['R_Ein'][i] * df['DEC_SN'][i] + (random.choice([1, -1])) * df['dec_off'][i]
    
    x[i] = xs
    y[i] = ys
    
    S.setLensPars(ml, rl, ql)
    S.setSourcePars(bs, ms, xs, ys, qs, ps, rs, 1)
    S.makeLens()

    try:
        magnification_gal[i] = S.magnification[1]
    except:
        magnification_gal[i] = 0

# Assign computed values back to DataFrame
df['x_gal'] = x
df['y_gal'] = y
df['Galaxy Magnification'] = magnification_gal

# Save updated DataFrame
df.to_csv('../outputs/SNIa/SNIa_GalaxyGalaxyPop_ideal.csv', index=False)
