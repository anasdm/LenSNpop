import numpy as np
import pandas as pd
from astropy.cosmology import Planck18
import sys
import scipy

lens_file = 'lens_catalogues_N100000.csv'
df_lenses = pd.read_csv('/users/anasdm/LenSNpop/' + lens_file)

C = 299792
zinterp = np.linspace(0, 2.5, 251)
distances = Planck18.angular_diameter_distance(zinterp).value
distfunc = scipy.interpolate.interp1d(zinterp, distances, fill_value='extrapolate', kind='cubic')

def R_Ein_precompute(sig, zl, z_source_range):
    ds = distfunc(z_source_range)
    dl = distfunc(zl)
    dls = (ds * (1 + z_source_range) - dl * (1 + zl)) / (1 + z_source_range)
    rein = np.power(sig, 2) * np.power((ds * np.power(C, 2)) / (4 * np.pi * dls), -1)
    rein[z_source_range <= zl] = 0  # Set to 0 where source redshift is less than lens redshift
    return rein * 206265

z_source_range = np.linspace(0, 2, 201)  # Adjust range and resolution as needed
einstein_radii_precomputed = np.array([R_Ein_precompute(sig, zl, z_source_range) for sig, zl in zip(df_lenses['sigma'], df_lenses['zl'])])

# Store precomputed values in the dataframe
df_lenses['einstein_radii_precomputed'] = list(einstein_radii_precomputed)



f_q = pd.read_csv('../data/SIE_SIS_ratio.csv')
interpolated_SIE_ratio = scipy.interpolate.interp1d(f_q['q'], f_q['Ratio'], fill_value = 'extrapolate', kind = 'cubic')
# Function to interpolate the Einstein radius
def interpolate_einstein_radius(precomputed_radii, z_source, z_source_range):
    interp_func = scipy.interpolate.interp1d(z_source_range, precomputed_radii, kind='cubic', fill_value='extrapolate')
    return interp_func(z_source)

def create_galaxy_galaxy_lenses(sourcefile, lensfile):
    df_sources = pd.read_csv('/users/anasdm/LenSNpop/' + sourcefile)

    lenses_per_source = 10
    lensed_df = pd.DataFrame()

    # Precompute some values for lenses
    lens_q = df_lenses['ql'].values
    lens_sig = df_lenses['sigma'].values
    lens_zl = df_lenses['zl'].values
    precomputed_radii = df_lenses['einstein_radii_precomputed'].values
    ratios = interpolated_SIE_ratio(lens_q)
    
    for source in df_sources.itertuples():  # Add tqdm for the progress bar
        # Einstein radius and weights calculations for all lenses for the current source
        source_z = source.z_host
        einstein_radii = np.array([interpolate_einstein_radius(precomputed_radii[i], source_z, z_source_range) for i in range(len(lens_zl))])
        lens_weights = ratios * einstein_radii ** 2

        # Select lenses based on weights
        indices = np.random.choice(len(df_lenses), lenses_per_source, p=lens_weights / lens_weights.sum())
        selected_lenses = df_lenses.iloc[indices].reset_index(drop=True)

        # Prepare source DataFrame for the selected lenses
        source_df = pd.DataFrame([source._asdict()] * lenses_per_source).reset_index(drop=True)

        # Combine lenses and sources
        combined_df = pd.concat([selected_lenses, source_df], axis=1, ignore_index=False)
        combined_df.reset_index(drop=True, inplace=True)

        # Compute additional required columns
        combined_df['rs'] = combined_df['R_d'] * Planck18.arcsec_per_kpc_comoving(combined_df['z_host'].values).value
        combined_df['R_Ein'] = einstein_radii[indices]

        # Append to the final DataFrame
        lensed_df = pd.concat([lensed_df, combined_df], ignore_index=True)

    return lensed_df

# Usage
SNtype = 'SNIa'
source_file = 'SNIa_source_catalogues_N10000.csv'

lens_galaxies_df=(create_galaxy_galaxy_lenses(source_file,lens_file))
lens_galaxies_df.to_csv('/users/anasdm/LenSNpop/outputs/%s_GalaxyGalaxyPop.csv'%SNtype)
