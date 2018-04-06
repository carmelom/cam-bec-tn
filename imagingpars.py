#!/usr/bin/python
#-*- coding: latin-1 -*-
"""Contains imaging parameters like effective pixelsize, absorption
coefficient, ..."""
class ImagingPars(object):
    """Base class for parameters of imaging system.

    @cvar description: descriptive name of settings, used for GUI selection.
    @cvar pixelsize: Size of area (in µm) which is imaged to one pixel of the cam.
    @cvar sigma0: cross section for light absorption
    """
    description = None
    pixelsize = 1
    sigma0 = 0
    expansion_time = 0
    mass = 0
    ODmax = 0 #maximum optical density

    def __str__(self):
        s = "%s, t_exp = %.1fms, OD_max = %.1f"%(
            self.description,
            self.expansion_time,
            self.ODmax)
        return s

class ImagingParsCMOS(ImagingPars):
    description = "cmos"
    pixelsize = 6.50/2.2#0.972  #pixelsize in µm x imaging magnification
    sigma0 = 1.5/3.14*(589e-9)**2
    mass = 23.0 * 1.66054e-27
    #palette = pylab.cm.gist_stern
    palette = "gist_stern"
        
class ImagingParsVerticalNa(ImagingPars):
    description = "vertical-Na"
    pixelsize = 4.40/0.5#0.972  #pixelsize in µm x imaging magnification
    sigma0 = 1.5/3.14*(589e-9)**2
    mass = 23.0 * 1.66054e-27
    #palette = pylab.cm.gist_stern
    palette = "gist_stern"

class ImagingParsAxialHRNa(ImagingPars):
    description = "axialHR-Na"
    pixelsize = 4.40/1.99 #0.972  #pixelsize in µm x imaging magnification
    sigma0 = 1.5/3.14*(589e-9)**2
    mass = 23.0 * 1.66054e-27
    #palette = pylab.cm.gist_stern
    palette = "gist_stern"

class ImagingParsAxialNa(ImagingPars):
    description = "axial-Na" #Magnif. last measurement 2018-04-03
    pixelsize = 4.40/1.053 #1.041 #0.972  #pixelsize in µm / imaging magnification
    sigma0 = (0.955) * 1.5/3.14*(589e-9)**2 # calibrated 2018-04-03
    mass = 23.0 * 1.66054e-27
    #palette = pylab.cm.gist_stern
    palette = "gist_stern"
    
class ImagingParsHorizontalHRNa(ImagingPars):
    description = "horizHR-Na"
    pixelsize = 4.40/2.042 #1.002  #pixelsize in µm x imaging magnification
    sigma0 = 1.5/3.14*(589e-9)**2
    mass = 23.0 * 1.66054e-27
    palette = "gist_stern"

class ImagingParsHorizontalNa(ImagingPars):
    description = "horiz-Na" #Magnif. last measurement 2018-04-03
    pixelsize = 4.40/1.3599 #1.362 #1.002  #pixelsize in µm / imaging magnification
    sigma0 = 1.5/3.14*(589e-9)**2
    mass = 23.0 * 1.66054e-27
    palette = "gist_stern"

class ImagingParsVerticalK(ImagingPars):
    description = "vert-K  "
    pixelsize = 4.40*2  #pixelsize in µm x imaging magnification
    sigma0 = 1.5/3.14*(767e-9)**2
    mass = 39.0 * 1.66054e-27
    palette = "spectral"

class ImagingParsHorizontalK(ImagingPars):
    description = "horiz-K  "
    pixelsize = 4.40*2  #pixelsize in µm x imaging magnification
    sigma0 = 1.5/3.14*(767e-9)**2
    mass = 39.0 * 1.66054e-27
    palette = "spectral"
