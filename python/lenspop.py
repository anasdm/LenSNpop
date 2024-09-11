#CLASS LENS POPULATION
#GENERATE REALISTIC POPULATION OF LENSES USING LENSPOP

import numpy as np
from astropy.units import Quantity
from astropy.table import Column
import pandas as pd
from scipy import interpolate
import pickle
import distances
import math
import scipy
from astropy.cosmology import Planck18


class RedshiftDependentRelation():
    def __init__(self,D=None,reset=True,cosmo=[0.31,0.69,Planck18.H(0).value]):
        self.beginRedshiftDependentRelation(D,reset=reset,cosmo=cosmo)

    def beginRedshiftDependentRelation(self,D,reset,zmax=10,cosmo=[0.3,0.7,0.7]):  
        self.zmax=zmax
        self.zbins,self.dz=np.linspace(0,self.zmax,401,retstep=True)
        self.z2bins,self.dz2=np.linspace(0,self.zmax,201,retstep=True)
        if D==None:
            D=distances.Distance(cosmo=cosmo)
        self.D=D
        
        if reset!=True:
            try:
            #load useful redshift splines
                splinedump=open("../data/LensPop/redshiftsplines.pkl","rb")
                self.Da_spline,self.Dmod_spline,self.volume_spline,self.Da_bispline=pickle.load(splinedump,encoding='latin1')
            except IOError or EOFError:   
                self.redshiftfunctions()
        else:
            self.redshiftfunctions()

    def redshiftfunctions(self):   
        D=self.D
        zbins=self.zbins
        z2bins=self.z2bins
        Dabins=zbins*0.0
        Dmodbins=zbins*0.0
        Da2bins=np.zeros((z2bins.size,z2bins.size))
        volumebins=zbins*0.0
        for i in range(zbins.size):
            Dabins[i]=D.Da(zbins[i])
            Dmodbins[i]=D.distance_modulus(zbins[i])
            volumebins[i]=D.volume(zbins[i])
        for i in range(z2bins.size):
            for j in range(z2bins.size):
                if j>i:
                    Da2bins[i,j]=D.Da(z2bins[i],z2bins[j])

        self.Da_spline=interpolate.splrep(zbins,Dabins)
        self.Dmod_spline=interpolate.splrep(zbins,Dmodbins)

        self.volume_spline=interpolate.splrep(zbins,volumebins)

        #z2d=iT.coords((z2bins.size,z2bins.size))*self.dz2
        #self.Da_bispline=interpolate.RectBivariateSpline(z2bins,z2bins,Da2bins)

        #pickle the splines
        #splinedump=open("../data/LensPop/redshiftsplines.pkl","wb")
        #pickle.dump([self.Da_spline,self.Dmod_spline,self.volume_spline,self.Da_bispline],splinedump,2)

    def Volume(self,z1,z2=None):
        if z2==None:
            return self.splev(z1,self.volume_spline)
        else:
            z1,z2=self.biassert(z1,z2)
            return self.splev(z2,self.volume_spline)-self.splev(z1,self.volume_spline)

    def Da(self,z1,z2=None,units="Mpc"):
        if units=="kpc":
            corfrac=1000
        elif units=="Mpc":
            corfrac=1
        else:
            print ("don't know those units yet")
        if z2==None:
            return self.splev(z1,self.Da_spline)*corfrac
        else:
            z1,z2=self.biassert(z1,z2)
            return self.Da_bispline.ev(z1,z2)*corfrac

    def Dmod(self,z):
        return self.splev(z,self.Dmod_spline)

    def splev(self,x,f_of_x_as_spline):
        return interpolate.splev(x,f_of_x_as_spline)

    def bisplev(self,x,y,f_ofxy_as_bispline):
        return interpolate.bisplev(x,y,f_ofxy_as_bispline)

    def biassert(self,z1,z2):
            try: len(z1)
            except TypeError:z1=[z1]
            try: len(z2)
            except TypeError:z2=[z2]
            if len(z1)==1 and len(z2)!=1:z1=np.ones(len(z2))*z1[0]
            if len(z2)==1 and len(z1)!=1:z2=np.ones(len(z1))*z2[0]
            assert len(z1)==len(z2)#,"get it together"
            return z1,z2

#====================================================================================


class EinsteinRadiusTools(RedshiftDependentRelation):
    def  __init__(self,D=None,reset=False):
        self.beginRedshiftDependentRelation(D,reset)
        self.c=299792

    def sie_sig(self,rein,zl,zs):
        self.c=299792
        ds=self.Da(zs)
        dls=self.Da(zl,zs)
        sig=(rein*(ds*self.c**2)/(206265*4*math.pi*dls))**0.5
        return sig
    def sie_rein(self,sig,zl,zs):
        self.c=299792
        ds=self.Da(zs)
        dls=self.Da(zl,zs)
        rein=sig**2*((ds*self.c**2)/(206265*4*math.pi*dls))**-1
        rein[rein<0]=0
        return rein


#====================================================================================
class Population(RedshiftDependentRelation):
    def  __init__(self):
        pass

    def draw_apparent_magnitude(self,M,z,band=None,colours=None):
        
        if band!=None:
            colours=self.colour(z,band)
        Dmods=self.Dmod(z)
        ml = M - colours + Dmods
        return ml
    
    def draw_apparent_size(self,r_phys,z): 
        rl = r_phys/(self.Da(z,units="kpc"))
        rl *= 206264 
        return rl
#======================
def E_eq(z):
        Om_m=Planck18.Om(0)
        Om_de=Planck18.Ode(0)
        return np.sqrt(Om_m*(1+z)**3+Om_de)

def dVol(z):
                if z==0:
                    return 0
                else:
                    c=3e5
                    D_H=c/Planck18.H(0).value
                    D=Planck18.angular_diameter_distance(z).value
                    return 4*np.pi*D_H*(1+z)**2*D**2/E_eq(z)
def phi(sigma):
            sigma[sigma==0]+=1e-6
            #phi_star=(8*10**-3)*self.D.h**3
            #alpha=2.32
            #beta=2.67
            #sigst=161 
            hr=(Planck18.H(0).value/(70))
            phi_star=(2.099*10**-2)*hr**3
            alpha= 0.94
            beta=1.85
            sigst=113.78
            phi=phi_star * \
                ((sigma*1./sigst)**alpha)*\
                np.exp(-(sigma*1./sigst)**beta)*beta/\
                math.gamma(alpha*1./beta)/\
                (1.*sigma)
            #phi*=(1+z)**(-2.5)
            return phi                
    
class Lenses(Population):
    def  __init__(self,Nlenses, zlmax=1.5,sig_min=100,D=None,reset=True,
                  bands=['F814W_ACS','g_SDSS','r_SDSS','i_SDSS','z_SDSS','Y_UKIRT','VIS']
                  ): 
            self.sig_min=sig_min
            self.zlmax=zlmax
            self.bands=bands
            self.N=Nlenses

            self.beginRedshiftDependentRelation(D,reset)

            with open('../data/LensPop/lenspopsplines.pkl', 'rb') as pickle_file:
                self.cdfdNdzasspline,self.cdfdNdsigz0asspline,self.cdfdsigdzasspline,self.dNdzspline,self.zlbins,self.zlmax,self.sigfloor,self.colourspline,self.bands= pickle.load(pickle_file, encoding='latin1')
    


            sig=np.linspace(100,400,301)
            zlbins,zl=np.linspace(0.01,1.5,201,retstep=True)
            sigbins=sig
            dNdzl=zlbins*0

            for i in range(len(zlbins)):
                z=zlbins[i]
                dphidsiggivenz=phi(sig)
                phisigspline=scipy.interpolate.splrep(sig,dphidsiggivenz)
                tot=scipy.interpolate.splint(100,400,phisigspline)
                if z!=0:
                    dNdzl[i]=tot*dVol(z)

            dNdzlf=scipy.interpolate.interp1d(zlbins,dNdzl,fill_value='extrapolate')

            Nofzcdf=np.cumsum(dNdzl)/np.sum(dNdzl)
            self.cdfdNdzasspline=scipy.interpolate.splrep(Nofzcdf,zlbins)

            dphidsig=phi(sig)
            cdfdNdsig=dphidsig.cumsum()/dphidsig.sum()
            self.cdfdNdsigspline=scipy.interpolate.splrep(cdfdNdsig,sig)

            #splinedump=open("../data/LensPop/lenspopsplines.pkl","rb")
            #self.cdfdNdzasspline,self.cdfdsigdzasspline,self.dNdzspline,self.zlbins,zlmax,sig_min,self.colourspline,bands=pickle.load(splinedump, encoding='latin1')
#Functions from LensPop.PopulationFunctions.py
    
    #def lensPopSplineDump(self):
    #        splinedump=open("../data/LensPop/lenspopsplines.pkl","wb")    	
    #        pickle.dump([self.cdfdNdzasspline,self.cdfdsigdzasspline,self.dNdzspline,self.zlbins,
    #        self.zlmax,self.sig_min,self.colourspline,self.bands],splinedump,2)
            
    def draw_sigma(self,N):
            sigs =interpolate.splev(np.random.random(N),self.cdfdNdsigspline)
            return sigs
    
    def draw_z(self,N):
            return interpolate.splev(np.random.random(N),self.cdfdNdzasspline)

    def draw_zsig(self,N):
            z=self.draw_z(N)
            sig=self.draw_sigma(N)
            return z,sig

    def EarlyTypeRelations(self,sigma,z=None,scatter=True,band=None):#z dependence not encoded currently
        #Hyde and Bernardi, M = r band absolute magnitude.
            V=np.log10(sigma)
            Mr=(-0.37+(0.37**2-(4*(0.006)*(2.97+V)))**0.5)/(2*0.006)
            if scatter:
                    Mr+=np.random.randn(len(Mr))*(0.15/2.4)

        #R=4.72+0.63*Mr+0.02*Mr**2 #rest-frame R_band size.
            R=2.46-2.79*V+0.84*V**2
            if scatter:
                    R+=np.random.randn(len(R))*0.11

        #convert to observed r band size;
            r_phys = 10**R

            return Mr,r_phys

    def colour(self,z,band):
            return interpolate.splev(z,self.colourspline[band])
            
    def Ndeflectors(self,z,zmin=0,fsky=1):
        if zmin>z:
            z,zmin=zmin,z
        N=interpolate.splint(zmin,z,self.dNdzspline)
        N*=fsky
        return N 

    def draw_flattening(self,sigma,z=None):
            x=sigma
            y=0.378-0.000572*x
            e=np.random.rayleigh(y)
            q=1-e
        #dont like ultraflattened masses:
            while len(q[q<0.2])>0 or len(q[q>1])>0:
                    q[q<0.2]=1-np.random.rayleigh(y[q<0.2])
                    q[q>1]=1-np.random.rayleigh(y[q>1])
            return q  #CHECK FOR np. INSTEAD OF NP. MAKE SURE IT IS ALL PYTHON3

    def drawLensPopulation(self):
            self.zl,self.sigl=self.draw_zsig(self.N)
            self.ql=self.draw_flattening(self.sigl)
            self.Mr,self.r_phys_nocol=self.EarlyTypeRelations(self.sigl,self.zl,scatter=True)
            self.ml={}
            self.rl={}
            self.r_phys={}
            for band in self.bands:
                    self.r_phys[band]=self.r_phys_nocol#could add a colorfunc here
                    if band !="VIS":
                        self.ml[band]=self.draw_apparent_magnitude(self.Mr,self.zl,band)
                    else: pass
                    self.rl[band]=self.draw_apparent_size(self.r_phys[band],self.zl) 
            return self.zl,self.sigl,self.ml,self.rl,self.ql
    

    def make_catalogue(self, save_df='yes'):
         df_lenses=pd.DataFrame(columns=['zl', 'sigma', 'ql', 'rl', 'mag_F814W_ACS', 'mag_g_SDSS','mag_r_SDSS', 'mag_i_SDSS', 'mag_z_SDSS', 'mag_Y_UKIRT'])
         df_lenses['zl']=self.zl
         df_lenses['sigma']=self.sigl
         df_lenses['ql']=self.ql
         df_lenses['rl']=self.rl['r_SDSS']
         #df_lenses['mag_F814W_ACS']=self.ml['F814W_ACS']
         df_lenses['mag_g_SDSS']=self.ml['g_SDSS']
         df_lenses['mag_r_SDSS']=self.ml['r_SDSS']
         df_lenses['mag_i_SDSS']=self.ml['i_SDSS']
         df_lenses['mag_z_SDSS']=self.ml['z_SDSS']
         #df_lenses['mag_Y_UKIRT']=self.ml['Y_UKIRT']

         if save_df=='yes':
              df_lenses.to_csv("../outputs/lens_catalogues_N%s.csv"%self.N,index=False)
		
         return df_lenses
#=========        	
if __name__=="__main__":
    L=Lenses(100000, zlmax=1.5)
    Lenses.drawLensPopulation(L)
    Lenses.make_catalogue(L)
