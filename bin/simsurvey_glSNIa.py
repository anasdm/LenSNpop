import simsurvey
import numpy as np
import pandas as pd
import json
from astropy.cosmology import Planck15
import sys
from scipy.interpolate import RectBivariateSpline
import os
#home_dir = os.getcwd()
import scipy
# Please enter the path to where you have placed the Schlegel, Finkbeiner & Davis (1998) dust map files
# You can also set the environment variable SFD_DIR to this path (in that case the variable below should be None)
sfd98_dir = '/users/anasdm/simsurvey/data/sfd98'
import sncosmo
#import simsurvey_tools as sst
#from scipy.interpolate import RectBivariateSpline as Spline2d
#from scipy.interpolate import UnivariateSpline as Spline1d
import pandas as pd
#from astropy import units as u
#from astropy import constants as const
#from scipy.interpolate import interp1d
#import extinction
from astropy.table import Table
import pickle
from extinction import fitzpatrick99

# %%
#SIMULATION

z_i,z_f=0.3,0.4 ###Simulated in \Delta 0.1 redshift bins to speed up simulation
ntransient=1000
dec_range=(-90,35)

def load_ztf_fields(filename='/users/anasdm/simsurvey/data/ZTF_Fields.txt', mwebv=False, galactic=False):
    """Load the ZTF fields propose by Eran from the included file.

    Parameters
    ----------
    filename: [str]
        File name of the ASCII file containing the field definitions

    mwebv: [bool]
        Include the Milky Way E(B-V) from the file in the output

    Return
    ------
    Dictionary of np.arrays with the field coordinates, IDs (and exinction)
    """
    fields = np.genfromtxt(filename, comments='%')

    out = {'field_id': np.array(fields[:,0], dtype=int),
           'ra': fields[:,1],
           'dec': fields[:,2]}

    if mwebv:
        out['mwebv'] = fields[:,3]
    if galactic:
        out['l'] = fields[:,4]
        out['b'] = fields[:,5]

    return out

def load_ztf_ccds(filename='/users/anasdm/simsurvey/data/ZTF_corners.txt', num_segs=16):
    """
    """
    ccd_corners = np.genfromtxt(filename, skip_header=1)/2
    ccds = [ccd_corners[4*k:4*k+4, :2] for k in range(num_segs)]
    return ccds

def load_ztf_filters():
    """
    """
    bands = {
        'ztfi' : 'data/ztfi_eff.txt',
        'ztfr' : 'data/ztfr_eff.txt',
        'ztfg' : 'data/ztfg_eff.txt',
    }

    for bandname in bands.keys() :
        fname = bands[bandname]
        b = np.loadtxt(fname)
        band = sncosmo.Bandpass(b[:,0], b[:,1], name=bandname)
        sncosmo.registry.register(band)



ccds=load_ztf_ccds()
dffields=pd.read_csv('/users/anasdm/simsurvey/data/lsstfields.csv')
fields={}
fields['field_id']=np.array(dffields['fieldID'],dtype=int)
fields['ra']=np.array(dffields['ra'])
fields['dec']=np.array(dffields['dec'])

df=pd.read_csv('/users/anasdm/simsurvey/surveyplans/lsstobs.csv')

obs = Table.from_pandas(df)

plan = simsurvey.SurveyPlan(time=obs['time'],
                            band=obs['band'],
                            obs_field=obs['field'],
                            skynoise=obs['skynoise'],
                            obs_ccd=None,
                            zp=obs['zp'],
                            comment=obs['comment'],
                            fields=fields,
                            ccds=ccds
                            )
mjd_range = (plan.pointings['time'].min(), plan.pointings['time'].max())

####DRAW STRETCH DISTRIBUTION##
meanx1=0.5
sigma_minus=1.3
sigma_plus=0.7
x1_range=np.linspace(-3,5)
def p_x1(x1):
    p_x=np.exp(-np.power(x1-meanx1,2)/(2*np.power(sigma_minus,2)))
    mask = x1>meanx1
    p_x[mask]=np.exp(-np.power(x1[mask]-meanx1,2)/(2*np.power(sigma_plus,2)))
    return p_x

cdf=np.cumsum(p_x1(x1_range))/np.sum(p_x1(x1_range))
cdfx1spline=scipy.interpolate.splrep(cdf, x1_range)

##AV DISTRIBUTION AS IN KESSLER ##
lambda_eff_z = 8922.78
def p_av(av):
    return 1/0.4 * np.exp(-av/0.4) + 1/(np.sqrt(2*np.pi)* 0.1)*np.exp(-np.power(av,2)/(0.01))

av_range = np.linspace(0,2)
cdf_av=np.cumsum(p_av(av_range))/np.sum(p_av(av_range))
cdfavspline=scipy.interpolate.splrep(cdf_av, av_range)

#INCLUDE STRETCH CORRELATIONS FROM GUY ET AL SALT2, 2007

dfSNe=pd.read_pickle('../outputs/SNIa/SNIa_larger_distribution.pkl')

p, w, f = sncosmo.read_griddata_fits(
    os.path.join(sncosmo.builtins.get_cache_dir(),
    'sncosmo/models/hsiao/Hsiao_SED_V3.fits')
)

class TimeSeriesSource_SL(sncosmo.Source):
    _param_names = ['amplitude', 's', 'mu1', 'mu2', 'mu3', 'mu4', 'dt1', 'dt2', 'dt3', 'index_cat']
    param_names_latex = ['A', 's', 'mu1', 'mu2', 'mu3', 'mu4', 'dt1', 'dt2', 'dt3', 'id']

    def __init__(self, phase, wave, flux, zero_before=True, zero_after=True, name=None,
                 version=None):
        self.name = name
        self.version = version
        self._zero_before = zero_before
        self._zero_after = zero_after
        self._phase = phase
        self._wave = wave
        self._parameters = np.array([1., 1., 1., 1., 1., 1., 0., 0., 0., 0])
        self._model_flux = RectBivariateSpline(self._phase, self._wave, flux, kx=3, ky=3)

    def _flux(self, phase, wave):
        f1 = self._model_flux(phase/ self._parameters[1], wave) * self._parameters[2]
        f2 = self._model_flux((phase-self._parameters[6])/ self._parameters[1], wave) * self._parameters[3]
        f3 = self._model_flux((phase-self._parameters[7])/ self._parameters[1], wave) * self._parameters[4]
        f4 = self._model_flux((phase-self._parameters[8])/ self._parameters[1], wave) * self._parameters[5]
        f = self._parameters[0] * (f1 + f2 + f3 + f4)

        if self._zero_before:
            mask = np.atleast_1d(phase) < self.minphase()
            f[mask, :] = 0.
        if self._zero_after:
            mask = np.atleast_1d(phase) > self.maxphase()
            f[mask, :] = 0.

        return f
        
source=TimeSeriesSource_SL(p,w,f)

dust = sncosmo.CCM89Dust()    

model = sncosmo.Model(source, effects=[dust, dust, dust], effect_names=['MW', 'lens', 'host'], effect_frames=['obs','free', 'rest'])

def random_parameters(redshifts, model, r_v=2., ebv_rate=0.11,
                          **kwargs):
    
    amplitude=[]
    mu1=[]
    mu2=[]
    mu3=[]
    mu4=[]

    t1=[]
    t2=[]
    t3=[]

    lensz=[]

    index=[]

    sampled_x1=scipy.interpolate.splev(np.random.random(len(redshifts)),cdfx1spline)
    stretch = np.abs(0.98 + (0.091 * sampled_x1) + (0.003 * np.power(sampled_x1, 2)) - (0.00075 * np.power(sampled_x1, 3)))
    
    sampled_av = scipy.interpolate.splev(np.random.random(len(redshifts)),cdfavspline)
    ebvhost = sampled_av/3.1

    for i, z in enumerate(redshifts):

        model2 = sncosmo.Model(source='hsiao')
        mask=np.abs(dfSNe['zs']-z)<=0.05

        df_glsn=dfSNe[mask]
        choose_glsn=df_glsn.sample(1)
        ind=choose_glsn.index[0]

        model2.set(z=z)
        
        t = np.linspace(-20,20, 41)
        model2.set_source_peakmag(float(choose_glsn['peakMag_z'].iloc[0]),'lsstz', 'ab') ##corrections to magnitude
        mag = model2.bandmag('lsstz', 'ab', t)
        amp = model2.get('amplitude')

        a_z = fitzpatrick99(np.array([lambda_eff_z]), sampled_av[i])[0]

        delta_mag = np.min(mag)-(float(choose_glsn['peakMag_z'].iloc[0])-a_z)
        factor = 10**(0.4*(delta_mag))
        amplitude.append(factor * amp)
  
        mu1.append(choose_glsn['mu1'][ind][0])
        mu2.append(choose_glsn['mu2'][ind][0])
        mu3.append(choose_glsn['mu3'][ind][0])
        mu4.append(choose_glsn['mu4'][ind][0])

        t1.append(choose_glsn['t2'][ind][0]/(1+z))
        t2.append(choose_glsn['t3'][ind][0]/(1+z))
        t3.append(choose_glsn['t4'][ind][0]/(1+z))

        lensz.append(float(choose_glsn['zl'].iloc[0]))
        index.append(ind)

    return {
            'amplitude': np.array(amplitude),
            's': np.array(stretch),
            'mu1': np.array(mu1),
            'mu2' : np.array(mu2),
            'mu3' : np.array(mu3),
            'mu4' : np.array(mu4),
            'dt1': np.array(t1),
            'dt2': np.array(t2),
            'dt3' : np.array(t3),
            'index_cat' : np.array(index),
            'lensz': np.array(lensz),
            'lensr_v': r_v * np.ones(len(redshifts)),
            'lensebv': np.random.exponential(ebv_rate, len(redshifts)),
            'hostr_v': 3.1 * np.ones(len(redshifts)),
            'hostebv': np.array(ebvhost)
            }

    
transientprop = dict(lcmodel=model,
                     lcsimul_func=random_parameters)

dec_range=(-90,30)
ra_range=(0,360)



tr = simsurvey.get_transient_generator((z_i, z_f),
                                      ntransient=1000,
                                      ratefunc=lambda z:3e-5,
                                      dec_range=dec_range,
                                      ra_range=ra_range,
                                      mjd_range=(mjd_range[0],
                                                 mjd_range[1]),
                                      transientprop=transientprop,
                                      sfd98_dir=sfd98_dir
                                      )

survey = simsurvey.SimulSurvey(generator=tr, plan=plan, n_det=3, threshold=5., sourcenoise=False)


lcs = survey.get_lightcurves(progress_bar=False, notebook=False)
lcs.save('../outputs/SNIa/simsurvey/simsurvey_ia.pkl')  
