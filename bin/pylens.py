import models
from math import pi
#import pylens
import MassModels

_PowerLawPars = [['b','eta','pa','q','x','y'],['b','eta','q','theta','x','y']]
_ExtShearPars = [['b','pa','x','y'],['b','theta','x','y']]
_MassSheetPars = [['b','x','y'],['b','x','y']]
_NFWPars = [['b','x','y'],['b','x','y']]

class PowerLaw(models._PowerLaw):
    """
    A subclass for power-law mass models. The `power-law' aspect doesn't
        currently work, but it does work for SIE models.
    """
    def __init__(self,name,var=None,const=None):
        if const is None:
            const = {}
        if var is None:
            var = {}
        # Check for all keys to be set
        keys = var.keys()+const.keys()
        keys.sort()
        if keys not in _PowerLawPars:
            import sys
            print ("Not all parameters defined!")
            sys.exit()
        models._PowerLaw.__init__(self)
        self.keys = keys
        self.values = {}
        self.vmap = {}
        for key in var.keys():
            self.values[key] = None
            self.vmap[key] = var[key]
        for key in const.keys():
            self.__setattr__(key,const[key])
        self.setPars()
        self.name = name
        self.NoFreeParams=False


    def __setattr__(self,key,value):
        if key=='pa':
            self.__dict__['pa'] = value
            if value is not None: 
                self.__dict__['theta'] = value*pi/180.
        elif key=='theta':
            if value is not None:
                self.__dict__['pa'] = value*180./pi
            self.__dict__['theta'] = value
        else:
            self.__dict__[key] = value

    def setPars(self):
        for key in self.vmap:
            self.__setattr__(key,self.vmap[key].value)


class SIE(PowerLaw):
    def __init__(self,name,var=None,const=None):
        c = {}
        for key in const.keys():
            c[key] = const[key]
        c['eta'] = 1.
        PowerLaw.__init__(self,name,var,c)


class ExtShear(models._ExtShear):
    def __init__(self,name,var=None,const=None):
        if const is None:
            const = {}
        if var is None:
            var = {}
        # Check for all keys to be set
        keys = var.keys()+const.keys()
        keys.sort()
        if keys not in _ExtShearPars:
            import sys
            print ("Not all parameters defined!",keys)
            sys.exit()
        models._ExtShear.__init__(self)
        self.keys = keys
        self.values = {}
        self.vmap = {}
        for key in var.keys():
            self.values[key] = None
            self.vmap[key] = var[key]
        for key in const.keys():
            self.__setattr__(key,const[key])
        self.setPars()
        self.name = name

    def __setattr__(self,key,value):
        if key=='pa':
            self.__dict__['pa'] = value
            if value is not None:
                self.__dict__['theta'] = value*pi/180.
        elif key=='theta':
            if value is not None:
                self.__dict__['pa'] = value*180./pi
            self.__dict__['theta'] = value
        else:
            self.__dict__[key] = value

    def setPars(self):
        for key in self.vmap:
            self.__setattr__(key,self.vmap[key].value)

class MassSheet(models._MassSheet):
    def __init__(self,name,var=None,const=None):
        if const is None:
            const = {}
        if var is None:
            var = {}
        # Check for all keys to be set
        keys = var.keys()+const.keys()
        keys.sort()
        if keys not in _MassSheetPars:
            import sys
            print("Not all parameters defined!",keys)
            sys.exit()
        models._MassSheet.__init__(self)
        self.keys = keys
        self.values = {}
        self.vmap = {}
        for key in var.keys():
            self.values[key] = None
            self.vmap[key] = var[key]
        for key in const.keys():
            self.__setattr__(key,const[key])
        self.setPars()
        self.name = name

    def __setattr__(self,key,value):
        self.__dict__[key] = value

    def setPars(self):
        for key in self.vmap:
            self.__setattr__(key,self.vmap[key].value)




class NFW(models._NFW):
    """
    A subclass for power-law mass models. The `power-law' aspect doesn't
        currently work, but it does work for SIE models.
    """
    def __init__(self,name,var=None,const=None):
        if const is None:
            const = {}
        if var is None:
            var = {}
        # Check for all keys to be set
        keys = var.keys()+const.keys()
        keys.sort()
        if keys not in _PowerLawPars:
            import sys
            print("Not all parameters defined!")
            sys.exit()
        models._PowerLaw.__init__(self)
        self.keys = keys
        self.values = {}
        self.vmap = {}
        for key in var.keys():
            self.values[key] = None
            self.vmap[key] = var[key]
        for key in const.keys():
            self.__setattr__(key,const[key])
        self.setPars()
        self.name = name

    def __setattr__(self,key,value):
        if key=='pa':
            self.__dict__['pa'] = value
            if value is not None: 
                self.__dict__['theta'] = value*pi/180.
        elif key=='theta':
            if value is not None:
                self.__dict__['pa'] = value*180./pi
            self.__dict__['theta'] = value
        else:
            self.__dict__[key] = value

    def setPars(self):
        for key in self.vmap:
            self.__setattr__(key,self.vmap[key].value)

def getDeflections(massmodels,points):
    if type(points)==type([]):
        x,y = points[0].copy(),points[1].copy()
    else:
        y,x = points[0].copy(),points[1].copy()
    if type(massmodels)!=type([]):
        massmodels = [massmodels]
    x0 = x.copy()
    y0 = y.copy()
    for massmodel in massmodels:
        xmap,ymap = massmodel.deflections(x,y)
        y0 -= ymap#.reshape(y.shape)#/scale
        x0 -= xmap#.reshape(x.shape)#/scale
    return x0.reshape(x.shape),y0.reshape(y.shape)



def lens_images(massmodels,sources,points,factor=1,getPix=False,csub=11):
    if type(points)==type([]):
        x,y = points[0].copy(),points[1].copy()
    else:
        y,x = points[0].copy(),points[1].copy()
    if type(massmodels)!=type([]):
        massmodels = [massmodels]
#    x0 = x.flatten()
#    y0 = y.flatten()
    x0 = x.copy()
    y0 = y.copy()
    for massmodel in massmodels:
        xmap,ymap = massmodel.deflections(x,y)
        y0 -= ymap#.reshape(y.shape)#/scale
        x0 -= xmap#.reshape(x.shape)#/scale
    x0,y0 = x0.reshape(x.shape),y0.reshape(y.shape)
   
    if type(sources)!=type([]):
        sources = [sources]
    out = x*0.
    for src in sources:
        out += src.pixeval(x0,y0,factor,csub=csub)

    if getPix==True:
        return out, x0,y0
    return out


def dblPlane(scales,massmodels,sources,points,factor,csub=11):
    if type(points)==type([]):
        x1,y1 = points[0].copy(),points[1].copy()
    else:
        y1,x1 = points[0].copy(),points[1].copy()
    out = x1*0.
    ax_1 = x1*0.
    ay_1 = x1*0.
    for l in massmodels[0]:
        xmap,ymap = l.deflections(x1,y1)
        ax_1 += xmap.reshape(ax_1.shape)
        ay_1 += ymap.reshape(ay_1.shape)
    for s in sources[0]:
        out += s.pixeval(x1-ax_1,y1-ay_1,factor,csub=csub)
    x2 = x1-scales[0,0]*ax_1
    y2 = y1-scales[0,0]*ay_1
    ax_2 = x2*0.
    ay_2 = y2*0.
    for l in massmodels[1]:
        xmap,ymap = l.deflections(x2,y2)
        ax_2 += xmap.reshape(ax_2.shape)
        ay_2 += ymap.reshape(ay_2.shape)
    for s in sources[1]:
        out += s.pixeval(x1-ax_1-ax_2,y1-ay_1-ay_2,factor,csub=csub)
    return out


def multiplePlanes(scales,massmodels,points):
    from numpy import zeros,eye,triu_indices
    if type(points)==type([]):
        x,y = points[0].copy(),points[1].copy()
    else:
        y,x = points[0].copy(),points[1].copy()
    out = x*0.
    nplanes = len(massmodels)
    tmp = scales.copy()
    scales = eye(nplanes)
    scales[triu_indices(nplanes,1)] = tmp
    ax = zeros((x.shape[0],x.shape[1],nplanes))
    ay = ax.copy()
    x0 = x.copy()
    y0 = y.copy()
    out = []
    for p in range(nplanes):
        for massmodel in massmodels[p]:
            xmap,ymap = massmodel.deflections(x0,y0)
            ax[:,:,p] += xmap.reshape(x.shape)
            ay[:,:,p] += ymap.reshape(y.shape)
        x0 = x-(ax[:,:,:p+1]*scales[:p+1,p]).sum(2)
        y0 = y-(ay[:,:,:p+1]*scales[:p+1,p]).sum(2)
        out.append([x0,y0])
    return out
