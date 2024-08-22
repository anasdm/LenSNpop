import numpy as np
from astropy.units import Quantity
from astropy.table import Column
import h5py
import pandas as pd

class Sources:

    def __init__(self, transient_type, Nsources):
        """
        This class defines the source population properties.

        :param transient_type: type of transient simulated. Options are ['AGN', 'KN', 'SLSN-I', 'SNII', 'SNIIb', 'SNIa', 'SNIb', 'SNIc', 'TDE'] -only Ia for now
        :param Nsources: total number of lens systems to be generated (int)
        ADD- :param scotch_dir: change directory for scotch catalogues
        """
        self.transient_type = transient_type
        self.Nsources = Nsources
        
    def generate_source_catalogue(self, zmin, zmax, save_df='yes'):
        """
        Generate SN+galaxies from the SCOTCH catalogues
        :param transient_type: options are 'Ia', 'II', 'Ib' and 'Ic'
        :param Nsources: total number of lens systems to be generated (int)
        :param z_min: minimum redshift for sources
        :param z_max: maximum redshift for sources
        :param save_df: option to save dataframe in outputs folder, defaults to 'yes'
        ADD- :param scotch_dir: change directory for scotch catalogue
        """
        try:
            scotch = h5py.File("../data/SourcePop/scotch_Z3_1.1.hdf5", "r")
        except OSError:
            print("Error: Unable to load SCOTCH catalog file.")

        transient_dict = {}
        hosts_dict = {}

        if self.transient_type == 'SNIa':
            scotch_transients = scotch['TransientTable']['SNIa']['SNIa-SALT2']
            scotch_hosts = scotch['HostTable']['SNIa']

        elif self.transient_type=='SNII':
            scotch_transients = scotch['TransientTable']['SNII']['SNII-Templates']
            scotch_hosts = scotch['HostTable']['SNII']

        elif self.transient_type=='SNIb':
            scotch_transients = scotch['TransientTable']['SNIb']['SNIb-Templates']
            scotch_hosts = scotch['HostTable']['SNIb']

        elif self.transient_type=='SNIc':
            scotch_transients = scotch['TransientTable']['SNIc']['SNIc-Templates']
            scotch_hosts = scotch['HostTable']['SNIc']

        elif self.transient_type=='SLSN':
            scotch_transients = scotch['TransientTable']['SLSN-I']['SLSN-I']
            scotch_hosts = scotch['HostTable']['SLSN-I']
            
        elif self.transient_type=='SNIIn':
            scotch_transients = scotch['TransientTable']['SNII']['SNIIn-MOSFIT']
            scotch_hosts = scotch['HostTable']['SNII']
        elif self.transient_type == 'KN':
            scotch_transients = scotch['TransientTable']['KN']['KN_K17']
            scotch_hosts = scotch['HostTable']['KN']
        elif self.transient_type == 'TDE':
            scotch_transients = scotch['TransientTable']['TDE']['TDE']
            scotch_hosts = scotch['HostTable']['TDE']
        else: 
            print('This SN type has not been included. The current options are "SNIa", "SNIb", "SNIc", "SNII" or "SLSN".')

            
            # Create a pandas dataframe for the transient properties
        for key in scotch_transients.keys():
            if len(np.shape(scotch_transients[key])) > 1:
                transient_dict[key] = [list(scotch_transients[key][i]) for i in range(len(scotch_transients[key][:50000+1, 0]))] #take fraction of schotch catalogue, 500000 first of total 2000000
            else:
                transient_dict[key] = scotch_transients[key][:50000+1]

            # Create a pandas dataframe for the host galaxy properties
        for key in scotch_hosts.keys():
            hosts_dict[key] = scotch_hosts[key][:]

            # Merge transient and host dictionaries into dataframes
        transient_df = pd.DataFrame(transient_dict)
        hosts_df = pd.DataFrame(hosts_dict) 
        merged_df = transient_df.merge(hosts_df, on=['TID', 'GID'], suffixes=('_transient', '_host'))

            # Quality cuts and filtering
        merged_df = merged_df[merged_df['logMstar'] < 900]
        merged_df = merged_df[merged_df['mag_i_host'] < 50]
        merged_df = merged_df[merged_df['z_host'] > zmin]
        merged_df = merged_df[merged_df['z_host'] < zmax]
        merged_df.reset_index(inplace=True, drop=True)
        print(len(merged_df))
            # Trim the merged dataframe to Nsources
        merged_df = merged_df.sample(self.Nsources)

        if save_df == 'yes':
            merged_df.to_csv("../outputs/%s_source_catalogues_N%s.csv" % (self.transient_type, self.Nsources), index=False)
 
        return merged_df

# Usage
if __name__ == "__main__":
    S = Sources('KN', 10000)  
    S.generate_source_catalogue(0.2, 2)

