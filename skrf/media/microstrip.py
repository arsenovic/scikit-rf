

'''
.. module:: skrf.media.microstrip
========================================
cpw (:mod:`skrf.media.microstrip`)
========================================

Coplanar waveguide class


This class was made from the technical documentation [#]_ provided
by the qucs project [#]_ , and from Mongia's book [#]_

.. [#] http://qucs.sourceforge.net/docs/technical.pdf
.. [#] http://www.qucs.sourceforge.net/
.. [#] "RF and Microwave Coupled-Line Circuits"
'''
from scipy.constants import  epsilon_0, mu_0
from scipy.special import ellipk
from numpy import real, imag,pi,sqrt,log,zeros, ones
from .media import Media
from ..tlineFunctions import skin_depth, surface_resistivity

class Microstrip(Media):
    '''
    Microstrip initializer

    A microstrip line can be created from either Z0 or w (strip width). 

    Parameters
    -------------
    frequency : :class:`~skrf.frequency.Frequency` object
        frequency band of the media
    Z0 : complex number, array-like, None
        desired characteristic impedance of microstrip in ohms. if None,
        `w` must be given
    w : number,  array-like, or None
        width of conductor, in m. unused if Z0 is given ,(it is 
        determined)
    ep_r : number, or array-like
        relative permativity of substrate
    h : number, or array-like
        thickness of substrate
    t : number, or array-like,
        conductor thickness, in m.
    rho: number, or array-like, optional
        resistivity of conductor (None)
    method: ['wheeler', 'hammerstad']
        method used to compute parameters
    z0 : number, array-like, or None
        the port impedance for media. Only needed if  its different
        from the characterisitc impedance of the transmission
    '''
    def __init__(self, frequency=None, Z0=None, w=None,ep_r=3, t=0, h=0,
                 rho=None, method='wheeler',z0=None, *args, **kwargs):
        Media.__init__(self, frequency=frequency,z0=z0)
        
        self.w = w
        self._Z0 = Z0
        self.ep_r = ep_r 
        self.h = h
        self.t = t
        
        self.rho = rho
                


    def __str__(self):
        f=self.frequency
        output =  \
                'Microstrip Waveguide Media.  %i-%i %s.  %i points'%\
                (f.f_scaled[0],f.f_scaled[-1],f.unit, f.npoints) + \
                '\n Z0 = %.2em, w= %.2em'% \
                (self.h, self.t)
        return output

    def __repr__(self):
        return self.__str__()

    @classmethod
    def Z0_from_w_static(self, w, ep_r,t,h):
        '''
        
        '''
        # eq 3.75 in Mongia
        dw = t* ((1+1/ep_r)/(2*pi) * \
                log( 10.872 / sqrt((t/h)**2 + ((1/pi)/(w/t+1.1))**2)))
        w_ = w + dw
        
        u = w_/h
        
        # eq 3.74 in Mongia
        Z0 = 42.4/(sqrt(ep_r+1)) * \
            log( 1 + (4/u) * ((14+8/ep_r)/11) * (4/u)+\
                 sqrt( ((14+8/ep_r)/11)**2 * (4/u)**2 + \
                       (1+1/ep_r)/2 * pi**2 ))
        
        return Z0

    @property
    def Z0(self):
        '''
        Characterisitc impedance
        '''
        Z0 =  self.Z0_from_w_static(w=self.w, ep_r=self.ep_r, t=self.t, h= self.h)
                       
        one = ones(len(self.frequency.f), dtype='complex')
        self._Z0 = Z0*one
        return self._Z0
    
    @property
    def epsilon_eff(self):
        Z0a = self.Z0_from_w_static(w=self.w, ep_r=1, t=self.t, h= self.h)
        Z0 =  self.Z0_from_w_static(w=self.w, ep_r=self.ep_r, t=self.t, h= self.h)
        
        return (Z0a/Z0)**2
        
    @property
    def gamma(self):
        '''
        Propagation constant
        
        '''
        k0 = self.frequency.w * sqrt(epsilon_0* mu_0)
        return 1j*sqrt(self.epsilon_eff*k0)
        
    def line(self,d,unit='deg',Z0=None, w=None,z0=None, embed = False, **kwargs):
        if (w is not None) and ( Z0 is not None):
            raise ValueError()
            
        if Z0 is not None:
            Z0_orig = self.Z0[:] 
            self.Z0 = Z0
            l = Media.line(self,d=d, unit=unit,z0=z0, embed=embed, **kwargs)
            self.Z0 = Z0_orig
        else:
            l = Media.line(self,d=d, unit=unit,z0=z0, embed=embed, **kwargs)
        
        
        if w is not None:
            w_orig = self.w
            self.w = w
            l = Media.line(self,d=d, unit=unit,z0=z0, embed=embed, **kwargs)
            self.w = w_orig
        else:
            l = Media.line(self,d=d, unit=unit,z0=z0, embed=embed, **kwargs)
            
        return l
