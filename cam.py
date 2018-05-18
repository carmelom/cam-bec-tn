#!/usr/bin/python2.7
#-*- coding: latin-1 -*-
"""Main file for camera program."""

from __future__ import with_statement




from sys import exit



import matplotlib
#matplotib.use('WxAgg')

import os, time, shutil, pickle, re, copy
import pylab
import numpy
import scipy.ndimage.filters
import scipy.signal
#import math

from numpy import ma

#gui
import wx, wx.aui, wx.grid
import wx.lib.delayedresult as delayedresult


from matplotlib.widgets import Button#, Cursor 
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wx import StatusBarWx

from matplotlib.backends.backend_wx import NavigationToolbar2Wx
#from toolbar import NavigationToolbar2Wx

from matplotlib.figure import Figure

#imports
from readsis import loadimg, loadimg3, read, write_raw_image
#loadimg = loadimg3

import settings
import fitting
import FitResultTableGrid
import DataPanel
import filewatch
import ding
from roi import ROI
import imagingpars
import ImageTree
import ImagePanel
from custom_events import *
#from profiling import Tic
import mean_sis


#force reload of modules (for development)
reload(FitResultTableGrid)
reload(settings)
reload(fitting)
reload(DataPanel)
reload(filewatch)
reload(ding)
reload(ImagePanel)
reload(ImageTree)


class CamCursor(object):
    """
    A horizontal and vertical line span the axes that and move with
    the pointer.  You can turn off the hline or vline spectively with
    the attributes

      horizOn =True|False: controls visibility of the horizontal line
      vertOn =True|False: controls visibility of the horizontal line

    And the visibility of the cursor itself with visible attribute
    """
    def __init__(self, ax, useblit=False, **lineprops):
        """
        Add a cursor to ax.  If useblit=True, use the backend
        dependent blitting features for faster updates (GTKAgg only
        now).  lineprops is a dictionary of line properties.  See
        examples/widgets/cursor.py.
        """
        self.ax = ax
        self.canvas = ax.figure.canvas

        self.canvas.mpl_connect('motion_notify_event', self.onmove)
        self.canvas.mpl_connect('draw_event', self.ondraw)

        self._visible = True
        self.horizOn = True
        self.vertOn = True
        self.useblit = useblit

        self.lineh = ax.axhline(ax.get_ybound()[0], visible=False, **lineprops)
        self.linev = ax.axvline(ax.get_xbound()[0], visible=False, **lineprops)

        self.background = None
        self.needclear = False
        
    def save_background(self):
        #print "save background"
        if self.useblit:
            self.background = self.canvas.copy_from_bbox(self.ax.bbox)

    def set_visible(self, visible):
        if visible == self._visible:
            return
        
        if visible == True:
            self.save_background()
        else:
            self.linev.set_visible(visible and self.vertOn)
            self.lineh.set_visible(visible and self.horizOn)

            pass
            #self._update()

        self._visible = visible
        

    def get_visible(self):
        return self._visible

    visible = property(get_visible, set_visible)

    def ondraw(self, event):
        #print "draw event"
        if self.visible:
            self.clear(event)
        else:
            pass
        #self.clear(event)

    def clear(self, event):
        'clear the cursor'
        #print "Cursor/clear"
        if self.useblit:
            self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        self.linev.set_visible(False)
        self.lineh.set_visible(False)

    def onmove(self, event):
        'on mouse motion draw the cursor if visible'
        if not self.visible: return #moved to beginning, this is the most important improvement

        if event.inaxes != self.ax:
            self.linev.set_visible(False)
            self.lineh.set_visible(False)

            if self.needclear:
                self._update() #instead of draw
                self.needclear = False
            return
        self.needclear = True
        self.linev.set_xdata((event.xdata, event.xdata))
        self.lineh.set_ydata((event.ydata, event.ydata))
        self.linev.set_visible(self.visible and self.vertOn)
        self.lineh.set_visible(self.visible and self.horizOn)

        self._update()


    def _update(self):
        if self.useblit:
            if self.background is not None:
                self.canvas.restore_region(self.background)
            else:
                #print "nothing to restore!"
                pass
            
            self.ax.draw_artist(self.linev)
            self.ax.draw_artist(self.lineh)
            self.canvas.blit(self.ax.bbox)
        else:

            self.canvas.draw_idle()

        return False



class Selector(object):
    "Base class for X/YSelector"
    
    def __init__(self, axis, num, tolerance=5, marker='+', color='y', alpha=0.5):
        self.axis = axis
        self.num = num
        self.tol = tolerance
        self.marker = marker
        self.color = color
        self.alpha = alpha

        self.cursor = CamCursor(self.axis, useblit=True)
        self.cursor.lineh.set_color(self.color)
        self.cursor.linev.set_color(self.color)
        self.init_cursor()
                             
        self.activated = False
        self.cursor.visible = False
        
    def init_cursor(self):
        pass
    
    def picker(self, artist, mousevent):
        for overaxis in [mousevent.inaxes]:
            if overaxis is self.axis:
                if not self.activated:
                    distance = self.distance(artist, mousevent)
                    if distance < self.tol: #hitted line
                        self.activated = True
                        #artist.set_alpha(1)
                        self.axis.figure.canvas.draw()
                        self.cursor.visible = True
                        return False, {}
                    
                else: #finish selection
                    self.activated = False
                    self.cursor.visible = False
                    newpos = self.update_position(artist, mousevent)
                    #artist.set_alpha(self.alpha)
                    overaxis.figure.canvas.draw() 
                    newpos['num'] = self.num
                    return True, newpos

        return False, {}

    def update_position(self, artist, mousevent):
        """set position of artist to position given by mouse and return
        new position as dict for use in picker callback"""
        pass

    def update_from_roi(self, roi):
        "set position from roi"
        pass

class YSelector(Selector):
    def __init__(self, axis, pos, num, **kwargs):
        Selector.__init__(self, axis, num, **kwargs)
        self.line = self.axis.axhline(pos,
                                      alpha=self.alpha,
                                      color=self.color,
                                      picker=self.picker)

    def init_cursor(self):
        self.cursor.vertOn = False

    def update_position(self, artist, mousevent):
        pos = round(mousevent.ydata)
        artist.set_ydata([pos] * 2)
        return dict(roiy=pos)

    def update_from_roi(self, roi):
        self.line.set_ydata([roi.gety(self.num)] * 2)
    
    def distance(self, artist, mousevent):
        linepos = artist.get_ydata()[0]
        #sx, sy     = self.axis.transData.xy_tup((0, linepos))
        sx, sy = self.axis.transData.transform((0, linepos))
        return abs(sy - mousevent.y)
            
        
class XSelector(Selector):
    def __init__(self, axis, pos, num, **kwargs):
        Selector.__init__(self, axis, num, **kwargs)
        self.line = self.axis.axvline(pos,
                                      alpha=self.alpha,
                                      color=self.color,
                                      picker=self.picker)
    def init_cursor(self):
        self.cursor.horizOn = False

    def update_position(self, artist, mousevent):
        pos = round(mousevent.xdata)
        artist.set_xdata([pos] * 2)
        return dict(roix=pos)

    def update_from_roi(self, roi):
        self.line.set_xdata([roi.getx(self.num)] * 2)

    def distance(self, artist, mousevent):
        linepos = artist.get_xdata()[0]
        #sx, sy   = self.axis.transData.xy_tup((linepos, 0))
        sx, sy = self.axis.transData.transform((linepos, 0))
        return abs(sx - mousevent.x)


class XYSelector(Selector):
    def __init__(self, axis, pos, num, **kwargs):
        Selector.__init__(self, axis, num, **kwargs)
        self.marker = self.axis.plot([pos[0]], [pos[1]],
                                     marker=self.marker,
                                     markersize=9,
                                     markeredgewidth=1.0,
                                     alpha=self.alpha,
                                     color=self.color,
                                     picker=self.picker)
        self.marker = self.marker[0]
        self.lineV = self.axis.axvline(pos[0],
                                        color=self.color,
                                        alpha=self.alpha)
        self.lineH = self.axis.axhline(pos[1],
                                        color=self.color,
                                        alpha=self.alpha)
        self.set_visible(False)
        
    def update_position(self, artist, mousevent):
        posx = round(mousevent.xdata)
        posy = round(mousevent.ydata)
        self.set_position((posx, posy))
        return dict(x=posx, y=posy)

    def distance(self, artist, mousevent):
        posx = artist.get_xdata()[0]
        posy = artist.get_ydata()[0]
        sx, sy = self.axis.transData.transform((posx, posy))
        return numpy.sqrt((sx - mousevent.x) ** 2 + (sy - mousevent.y) ** 2)

    def set_position(self, point):
        posx, posy = point
        self.marker.set_xdata([posx])
        self.marker.set_ydata([posy])
        self.lineV.set_xdata([posx] * 2)
        self.lineH.set_ydata([posy] * 2)

    def get_position(self):
        posx = self.marker.get_xdata()[0]
        posy = self.marker.get_ydata()[0]
        return (posx, posy)

    position = property(get_position, set_position)

    def set_visible(self, visible):
        self.cross_visible = visible
        self.lineV.set_visible(self.cross_visible)
        self.lineH.set_visible(self.cross_visible)

    def get_visible(self):
        visible = self.lineV._visible
        return visible


class CamSplashScreen(wx.SplashScreen):
    def __init__(self):
        bitmap = wx.Image(os.path.join(settings.bitmappath, 'cam_splash.png'),
                          wx.BITMAP_TYPE_PNG).ConvertToBitmap()
        wx.SplashScreen.__init__(self,
                                 bitmap,
                                 wx.SPLASH_CENTRE_ON_SCREEN | wx.SPLASH_TIMEOUT,
                                 1000, None, - 1)

class ImgPanel(wx.Panel):
    """Class for displaying images and handling fit results, realized
    as wxPanel"""

    ID_Autoscale = wx.NewId()
    
    def __init__(self, parent, panel_id, img, region_of_interest, fitroutine, palette=pylab.cm.jet):

        #Initialize wx stuff 
        wx.Panel.__init__(self,
                          parent,
                          size=(550, 350))

        self.fig = pylab.Figure()
        self.canvas = FigureCanvas(self, -1, self.fig)
        
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.SetSizer(self.sizer)
        self.Fit()

        #setup (and modify) toolbar
        self.toolbar = NavigationToolbar2Wx(self.canvas)
        self.toolbar.AddLabelTool(self.ID_Autoscale,
                                  'Autoscale',
                                  wx.Bitmap(os.path.join(settings.bitmappath,
                                                         'autoscale.png'),
                                            wx.BITMAP_TYPE_PNG),
                                  shortHelp = 'Autoscale',
                                  longHelp = 'automatic scaling')
        wx.EVT_TOOL(self, self.ID_Autoscale, self.OnAutoscale)
        self.toolbar.Realize()
        tw, th = self.toolbar.GetSizeTuple()
        fw, fh = self.canvas.GetSizeTuple()
        self.toolbar.SetSize(wx.Size(fw, th))
        self.sizer.Add(self.toolbar, 0, wx.LEFT | wx.EXPAND)
        self.toolbar.update()

        #store input data
        self.id = panel_id
        self.rawimg = img
        self.img = img #for the moment
        self.fit = fitroutine
        self.roi = region_of_interest
        
        #initialize attributes
        self.x = range(img.shape[0])
        self.y = range(img.shape[1])
        self.fitpars = fitting.FitParsNoFit()
        self.compensate_offset = False
        self.target_results = None
        self.difference=False
        self.filter=False
        self.filter_scalar=1.0
        self.filter_sigma=1.0
        
        self.reference_diff=False
        self.reference_img=None
        
        self.rotation_angle=0.0
        self.rotation_active=False
        self.bec_residuals=None
        
        #create axes
        #self.axtop   = self.fig.add_axes([0.2, 0.925, 0.6, 0.05]) #top (reload)
        self.aximg = self.fig.add_axes([0.2, 0.3, 0.6, 0.65])    #image
        self.axvprof = self.fig.add_axes([0.05, 0.3, 0.13, 0.65], sharey=self.aximg) #vertical profile
        self.axhprof = self.fig.add_axes([0.2, 0.05, 0.6, 0.23], sharex=self.aximg) #horizontal profile
        #self.axclrb = self.fig.add_axes([0.82, 0.3, 0.05, 0.65])  #colorbar
        #self.axtxt = self.fig.add_axes([0.82, 0.05, 0.15, 0.23])#text
        self.axclrb = self.fig.add_axes([0.82, 0.05, 0.05, 0.23]) #colorbar
        self.axtxt = self.fig.add_axes([0.82, 0.3, 0.05, 0.65])  #text


        #text area
        self.axtxt.set_axis_off()
        self.htxt = self.axtxt.text(0, 0, "Hi",
                                    axes=self.axtxt,
                                    multialignment='left',
                                    transform=self.axtxt.transAxes,
                                    family='monospace',
                                    )
        
        #initialize image area
        
        self.himg = self.aximg.imshow(img, cmap=palette, vmin= - 0.05, vmax=1.5,
                                      origin="upper",
                                      interpolation='bicubic' #'nearest', 'spline36'
                                      )
        
        self.aximg.set_axis_off()
        self.hclrbr = pylab.colorbar(self.himg, self.axclrb)


        #initialize profiles
        self.create_profiles()

        #initialize contour handles
        self.hcontours = []

        #create borders ROI
        self.hrx = []

        self.hrx.append(
            XSelector(
                self.aximg,
                self.roi.xmin,
                num=0))

        self.hrx.append(
            XSelector(
                self.aximg,
                self.roi.xmax,
                num=1))
        
        self.hry = []
        self.hry.append(
            YSelector(
                self.aximg,
                self.roi.ymin,
                num=0
            ))

        self.hry.append(
            YSelector(
                self.aximg,
                self.roi.ymax,
                num=1
            ))
        
        #initialize markers
        self.markers = []
        self.markers.append(
            XYSelector(self.aximg, (440, 295), num=0, marker='x',color='w', alpha=0.8))
        self.markers.append(
            XYSelector(self.aximg, (485, 295), num=1, color='w', alpha=0.8))
        self.markers.append(
            XYSelector(self.aximg, (440, 340), num=2, color='g', alpha=0.8))


        #connect event handler for manipulation of ROI-borders
        self.fig.canvas.mpl_connect('pick_event', self.onpick)

        #set initialize limits for image _and_ profiles
        self.center_roi()

        #TODO: queak
        self.toolbar.push_current()

        self.Bind(wx.EVT_ERASE_BACKGROUND, self.OnEraseBackground)
        
        self.fitjobID = 0
        self.Bind(EVT_FIT_COMPLETED, self.OnFitCompleted)
        self.update()

        self.axhprof.autoscale_view(scalex = False, scaley = True)
        self.axvprof.autoscale_view(scalex = True, scaley = False)

    def OnEraseBackground(self, event):
        event.Skip()
        
    def create_profiles(self):
        r = self.roi
        imgroi = self.img[r.y, r.x]
        self.hprof = imgroi.sum(0)
        self.vprof = imgroi.sum(1)

        #self.hprof = ma.array(self.hprof, mask = ~numpy.isfinite(self.hprof))
        #self.vprof = ma.array(self.vprof, mask = ~numpy.isfinite(self.vprof))
        #takex = numpy.isfinite(self.hprof)
        #takey = numpy.isfinite(self.vprof)

        x = r.xrange_clipped(self.img)
        y = r.yrange_clipped(self.img)

        self.hhprof = self.axhprof.plot(x, self.hprof, zorder=10,
                                        #drawstyle = 'steps-mid'
                                        )[0]
        self.hvprof = self.axvprof.plot(self.vprof, y, zorder=10,
                                        #drawstyle = 'steps-mid'
                                        )[0]
        
        self.hhproffit = []
        self.hvproffit = []
        for color, zorder in zip(('r', 'm', 'k:'), (5, 4, 0)):
            self.hhproffit += self.axhprof.plot([], [], color, zorder=zorder)
            self.hvproffit += self.axvprof.plot([], [], color, zorder=zorder)

    def update_profiles(self):
        r = self.roi
        imgroi = self.img[r.y, r.x]

        if self.compensate_offset:
            cimg = imgroi - self.img_background
        else:
            cimg = imgroi

        self.hprof = cimg.filled().sum(0)
        self.vprof = cimg.filled().sum(1)
        
        self.hhprof.set_data(r.xrange_clipped(self.img), self.hprof)
        self.hvprof.set_data(self.vprof, r.yrange_clipped(self.img))
        
        
        
        
        
        
    def mpl_imsave(self,fname, arr, vmin=None, vmax=None, cmap=None, format=None, origin=None, dpi=100):
        """
        Save an array as in image file.
        
        The output formats available depend on the backend being used.
        
        Arguments:
        *fname*:
        A string containing a path to a filename, or a Python file-like object.
        If *format* is *None* and *fname* is a string, the output
        format is deduced from the extension of the filename.
        *arr*:
        An MxN (luminance), MxNx3 (RGB) or MxNx4 (RGBA) array.
        Keyword arguments:
        *vmin*/*vmax*: [ None | scalar ]
        *vmin* and *vmax* set the color scaling for the image by fixing the
        values that map to the colormap color limits. If either *vmin*
        or *vmax* is None, that limit is determined from the *arr*
        min/max value.
        *cmap*:
        cmap is a colors.Colormap instance, eg cm.jet.
        If None, default to the rc image.cmap value.
        *format*:
        One of the file extensions supported by the active
        backend. Most backends support png, pdf, ps, eps and svg.
        *origin*
        [ 'upper' | 'lower' ] Indicates where the [0,0] index of
        the array is in the upper left or lower left corner of
        the axes. Defaults to the rc image.origin value.
        *dpi*
        The DPI to store in the metadata of the file. This does not affect the
        resolution of the output image.
        """
        from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
        from matplotlib.figure import Figure
    
        figsize = [x / float(dpi) for x in (arr.shape[1], arr.shape[0])]
        fig = Figure(figsize=figsize, dpi=dpi, frameon=False)
        canvas = FigureCanvas(fig)
        im = fig.figimage(arr, cmap=cmap, vmin=vmin, vmax=vmax, origin=origin)
        fig.savefig(fname, dpi=dpi, format=format, transparent=True)
    
    
    def create_bitmap(self):
        if self.difference and self.imgfit_tot!= None:
            outimg=self.img-self.imgfit_tot
        elif self.compensate_offset:
            outimg=self.img - self.img_background
        else:
            outimg=self.img
        
        
        if self.reference_diff and (self.reference_img!=None):
            outimg-=self.reference_img
        if self.filter:
            outimg=scipy.ndimage.filters.gaussian_filter(outimg,[float(self.filter_sigma),float(self.filter_sigma)])
            #outimg=scipy.signal.wiener(outimg,[float(self.filter_sigma),float(self.filter_sigma)])
            #outimg=scipy.ndimage.filters.median_filter(outimg,[float(self.filter_sigma),float(self.filter_sigma)])
        outimg*=self.filter_scalar
        
        self.calc_residuals()
        
        return outimg
    
    
    def calc_residuals(self):
        if isinstance(self.fitpars, (fitting.FitParsBimodal2d, fitting.FitParsTF2d)):
            
            if self.reference_diff and self.reference_img!=None:
                myimgroi = self.img - self.reference_img
            else:
                myimgroi=self.img
                
            if self.filter:
                myimgroi=scipy.ndimage.filters.gaussian_filter(myimgroi,[float(self.filter_sigma),float(self.filter_sigma)])
            
            r = self.roi
            myimgroi=myimgroi[r.y,r.x]
            
            if isinstance(self.fitpars, fitting.FitParsBimodal2d):
                #imgfitgauss = (self.imgfit[1] - self.img_background)
                imgfitbec   = (self.imgfit[0] - self.imgfit[1])
                mymask = numpy.ma.getmask(numpy.ma.masked_less_equal(imgfitbec, numpy.ma.maximum(imgfitbec)*0.00001))
                myimgfitbec = numpy.ma.masked_array(data=self.imgfit[0], mask=mymask)
            elif isinstance(self.fitpars, fitting.FitParsTF2d):
                imgfitbec = self.imgfit[0]  - self.img_background
                myimgfitbec = self.imgfit[0]
                
            max = numpy.ma.maximum(myimgfitbec)
            min = numpy.ma.minimum(myimgfitbec)
            
            masked_fit = numpy.ma.masked_less_equal(myimgfitbec, (max-min)*0.2+min)
            numpy.ma.set_fill_value(masked_fit, numpy.nan)
            numpy.ma.set_fill_value(myimgroi, numpy.nan)
            self.bec_residuals = numpy.ma.mean(((myimgroi - masked_fit)/imgfitbec)**2)
        else:
            self.bec_residuals=None
        
        self.update_parameters()
        
        return self.bec_residuals
    
    
    
    def save_bitmap(self,filename,cmap):
        min,max=self.himg.get_clim()
        self.mpl_imsave(filename, self.create_bitmap(), cmap=cmap,vmin= min, vmax=max)

    def clear_image_axis_obsolete(self):
        """clean image axis from old contours. for this, clear axis
        completely and readd artist (cursors, markers, ...)"""
       
        self.aximg.clear()
        self.aximg.set_axis_off()
        self.aximg.add_artist(self.himg)
        for marker in self.hrx + self.hry:
            self.aximg.add_artist(marker.cursor.lineh)
            self.aximg.add_artist(marker.cursor.linev)
            self.aximg.add_artist(marker.line)
        for marker in self.markers:
            self.aximg.add_artist(marker.marker)
        
        self.aximg.set_autoscale_on(False)

    def clear_contours(self):
        """remove contour lines from image axis"""
        self.aximg.collections = [h for h in self.aximg.collections if h not in self.hcontours]
        self.hcontours = []

    def update_contours(self):
        try:
            self.clear_contours()
            if self.show_contours:
                xlim = numpy.copy(self.aximg.get_xlim())
                ylim = numpy.copy(self.aximg.get_ylim())

                if isinstance(self.fitpars, fitting.FitParsBimodal2d):
                    imgfitgauss = (self.imgfit[1] - self.img_background)
                    imgfitbec   = (self.imgfit[0] - self.imgfit[1])
                    cs = self.aximg.contour(self.roi.xrange_clipped(self.img),
                               self.roi.yrange_clipped(self.img),
                               imgfitgauss,
                               [numpy.exp(-2)*imgfitgauss.max(),],
                               colors = 'w',
                               alpha = 0.5,
                               linewidths = 1.0,
                               linestyle = 'dashed')

                    cs2 = self.aximg.contour(self.roi.xrange_clipped(self.img),
                               self.roi.yrange_clipped(self.img),
                               imgfitbec,
                               [0.05*imgfitbec.max(),],
                               colors = 'w',
                               alpha = 0.5,
                               linewidths = 1.0)#
                    cs2.collections[0].set_linestyle('dashed')
                    self.hcontours = cs.collections + cs2.collections
                    
                elif isinstance(self.fitpars, fitting.FitParsBimodalIntegrate1dSplit):
                    self.hcontours = []

                elif isinstance(self.fitpars, fitting.FitParsGauss2d):
                    cs = self.aximg.contour(self.roi.xrange_clipped(self.img),
                               self.roi.yrange_clipped(self.img),
                               self.imgfit[0],
                               [0.37*self.fitpars.OD],
                               colors = 'w',
                               alpha = 0.5,
                               linewidths = 1.0)#
                    self.hcontours = cs.collections

                self.aximg.set_xlim(xlim)
                self.aximg.set_ylim(ylim)

            else:
                self.hcontours = []
                
        except Exception, e:
            print "Error in update contours: ", e
            import traceback
            traceback.print_exc()

    def update(self):
        "calculate profiles and fit results, perform redraw"
        self.invalidate_profiles_fit()
        self.invalidate_parameters()
        self.bec_residuals = None
        self.clear_contours()
        self.create_image()
        self.himg.set_data(self.img)
        self.redraw()
        
        self.fit.rotation_active=self.rotation_active
        self.fit.rotation_angle=self.rotation_angle

        self.fitjobID += 1
        delayedresult.startWorker(consumer = self,
                                  workerFn = self._fit_result_producer,
                                  wargs = (self.fit, #worker arguments
                                           self.img.copy(),
                                           copy.copy(self.roi)), 
                                  cargs = (FitCompletedEvent,), #consumer arguments
                                  ckwargs = {'resultAttr': 'fitresults',
                                             'target': self.target_results}, 
                                  jobID = self.fitjobID)
        
    def _fit_result_producer(self, fit, img, roi):
        if self.difference:
            return fit.do_fit(img, roi, True)
        else:
            return fit.do_fit(img, roi)

    def OnFitCompleted(self, event):
        fitresult = event.fitresults
        jobID = fitresult.getJobID()
        if jobID != self.fitjobID:
            print "Received fit results from outdated fit. Ignoring these results."
            return
        try:
            self.imgfit, self.img_background, self.fitpars, self.imgfit_tot = fitresult.get()
        except Exception, e:
            #if fit failed, create save settings, don't update profiles and
            #contours
            print "Error in multi-threaded fit:", e
            self.imgfit = numpy.empty(shape=(0,0))
            self.img_background = numpy.array([0])
            self.fitpars = fitting.FitParsNoFit()
        else:
            self.update_profiles_fit()
            self.update_contours()

        self.update_profiles()
        self.update_parameters()
        self.update_image()
        self.redraw()

        #tell the world that fit is completed
        if jobID > 1:
            evt = FitResultsEvent(id=0,
                                  source = self.id,
                                  target = event.target,
                                  fitpars = self.fitpars)
            wx.PostEvent(self.canvas, evt)

    def analyze_image(self, img, target):
        """show image and perform fit. store target information for
        later use."""
        self.rawimg = img
        self.target_results = target
        self.update()
        
    def create_image(self):
        """
        apply image filtering, create self.img, (don't update image display data)
        """
        img = self.rawimg
        
        ##compensate for finite optical density
        ODmax = self.fit.imaging_pars.ODmax #TODO: this might fail if
                                            #fit class is changed.
        if ODmax > 0:
            img = numpy.log((1 - numpy.exp(- ODmax)) / (numpy.exp(- img) - numpy.exp(- ODmax)))
            #TODO: remove invalid entries
            #TODO: this should not happen here!
        
        self.img = ma.array(img, mask= ~ numpy.isfinite(img))
        if ODmax > 0:
            self.img.set_fill_value(ODmax)
        else:
            self.img.set_fill_value(3) #which to take?
        
    def redraw(self):   
        self.fig.canvas.draw()
        
    def invalidate_parameters(self):
        self.htxt.set_alpha(0.6)
        self.htxt.set_color('k')
        self.htxt.set_text('fitting...')

    def update_parameters(self):
        s = unicode(self.fitpars)
        if self.bec_residuals!=None:
            s+="\nBEC residuals:\n%.5f"%self.bec_residuals
            
        self.htxt.set_text(s)

        if self.fitpars.valid:
            self.htxt.set_color('k')
            self.htxt.set_alpha(1)
        else:
            self.htxt.set_color('r')
            self.htxt.set_alpha(0.2)
        
    def invalidate_profiles_fit(self):
        for h in self.hhproffit + self.hvproffit:
            h.set_data([], [])

    def update_profiles_fit(self):
        for k, img in enumerate(self.imgfit):
            if self.compensate_offset:
                cimg = img - self.img_background
            else:
                cimg = img

            self.hhproffit[k].set_data(self.roi.xrange_clipped(self.img), cimg.sum(0))
            self.hvproffit[k].set_data(cimg.sum(1), self.roi.yrange_clipped(self.img))

        if self.compensate_offset:
            self.hhproffit[ - 1].set_data([0, 1250], [0, 0])
            self.hvproffit[ - 1].set_data([0, 0], [0, 1040])

        
        
        
    #def do_fit(self):
    #    self.imgfit, self.img_background, self.fitpars = self.fit.do_fit(self.img, self.roi)

    def update_image(self):
        self.himg.set_data(self.create_bitmap())

                
    def set_compensate_offset(self, value=True):
        self.compensate_offset = value
        self.update() #TODO: some optimazation?
        #self.update_profiles()
        #self.invalidate_profiles_fit()
        #self.update_profiles_fit()
        #self.update_image()
        #self.redraw()

    def set_marker_guess(self, value=False):
        self.marker_guess = value
        self.update() 

    def set_roi(self, roi):
        "set roi, update selectors to match roi"
        self.roi = roi

        for selector in self.hrx + self.hry:
            selector.update_from_roi(roi)

        #TODO: update fits?

    def center_roi(self, border=50):
        self.axhprof.set_xlim(self.roi.xmin - border, self.roi.xmax + border)
        self.axvprof.set_ylim(self.roi.ymax + border, self.roi.ymin - border)

    def onpick(self, event):
        need_update = False
        
        "handle pick events. Knows about handling events generated by ROI-picks"
        try:
            self.roi.setx(event.roix, event.num)
        except AttributeError:
            pass
        else:
            need_update = True

        try:
            self.roi.sety(event.roiy, event.num)
        except AttributeError:
            pass
        else:
            need_update = True

        try:
            x, y = event.x, event.y
        except AttributeError:
            pass
        else:
            self.redraw()
        
        if need_update:
            self.update()

    def OnAutoscale(self, event):
        prof = self.vprof
        mi = prof.min()
        ma = prof.max()
        if self.compensate_offset:
            mi = 0
        d = 0.2*(ma - mi)
        self.axvprof.set_xlim(mi-d, ma+d)
        
        prof = self.hprof
        mi = prof.min()
        ma = prof.max()
        if self.compensate_offset:
            mi = 0
        d = 0.2*(ma - mi)
        self.axhprof.set_ylim(mi-d, ma+d)
        
        border = 50
        self.axhprof.set_xlim(self.roi.xmin - border, self.roi.xmax + border)
        self.axvprof.set_ylim(self.roi.ymax + border, self.roi.ymin - border)
        
        self.redraw()
    
##results name

class ImgAppAui(wx.App):

    #Menu: View
    ID_ShowNa = wx.NewId()
    ID_ShowK = wx.NewId()
    ID_ShowTools = wx.NewId()
    ID_ShowSavedImages = wx.NewId()
    
    #Menu: Perspectives
    ID_Perspective = wx.NewId()
    ID_PerspectiveFull = wx.NewId()
    ID_PerspectiveNa = wx.NewId()
    ID_PerspectiveK = wx.NewId()
    
    #Menu: Fit
    ID_Reload = wx.NewId()
    
    ID_FitShowContours = wx.NewId()
    
    ID_FitNaNone = wx.NewId()
    ID_NaImageIntegration = wx.NewId()
    ID_FitNaGauss1D = wx.NewId()
    ID_FitNaLorentzGauss1D = wx.NewId()
    ID_FitNaGaussGauss1D = wx.NewId()
    ID_FitNaGauss = wx.NewId()
    ID_FitNaGaussSym = wx.NewId()
    ID_FitNaGaussBose = wx.NewId()
    ID_FitNaTF = wx.NewId()
    ID_FitNaTFSlice = wx.NewId()
    ID_FitNaBimodal = wx.NewId()
    ID_FitNaBimodalSplit = wx.NewId()
    ID_FitNaThomasFermiHole = wx.NewId()
    ID_FitNaBimodalIntegrate1dSplit = wx.NewId()
    ID_FitNaGaussIntegrate1dPinning = wx.NewId()
    ID_FitNaBoseBimodal = wx.NewId()
    ID_FitNaBimodalGaussGauss = wx.NewId()
    ID_FitNaNGauss = wx.NewId()
    
    ID_FitKNone = wx.NewId()
    ID_KImageIntegration = wx.NewId()
    ID_FitKGauss1D = wx.NewId()
    ID_FitKLorentzGauss1D = wx.NewId()
    ID_FitKGaussGauss1D = wx.NewId()
    ID_FitKGauss = wx.NewId()
    ID_FitKGaussSym = wx.NewId()
    ID_FitKGaussBose = wx.NewId()
    ID_FitKTF = wx.NewId()
    ID_FitKTFSlice = wx.NewId()
    ID_FitKBimodal = wx.NewId()
    ID_FitKBimodalSplit = wx.NewId()
    ID_FitKThomasFermiHole = wx.NewId()
    ID_FitKBimodalIntegrate1dSplit = wx.NewId()
    ID_FitKGaussIntegrate1dPinning = wx.NewId()
    ID_FitKBimodalGaussGauss = wx.NewId()
    ID_FitKBoseBimodal = wx.NewId()

    
    ID_CompensateOffsetNa = wx.NewId()
    ID_CompensateOffsetK = wx.NewId()
    
    ID_MarkerGuessNa = wx.NewId()
    ID_MarkerGuessK = wx.NewId()

    #Menu: Display settings
    ID_K_ContrastHigh = wx.NewId()
    ID_K_ContrastUltra = wx.NewId()
    ID_K_ContrastNormal = wx.NewId()
    ID_K_ContrastLow = wx.NewId()
    
    ID_Na_ContrastHigh = wx.NewId()
    ID_Na_ContrastUltra = wx.NewId()
    ID_Na_ContrastNormal = wx.NewId()
    ID_Na_ContrastLow = wx.NewId()
    
    ID_MarkerA = wx.NewId()
    ID_MarkerB = wx.NewId()
    ID_MarkerC = wx.NewId()
    
    #Menu: Results
    ID_ResultsPlot = wx.NewId()
    ID_ResultsSave = wx.NewId()
    ID_ResultsSaveAs = wx.NewId()
    ID_ResultsLoad = wx.NewId()
    ID_ResultsSaveTemplate = wx.NewId()
    ID_ResultsApplyTemplate = wx.NewId()
    ID_ResultsNew = wx.NewId()
    
    #Menu: Settings
    ID_SettingsSave = wx.NewId()
    ID_SettingsLoad = wx.NewId()
    
    #Menu: Autosav
    ID_AutoSaveAbs = wx.NewId()
    ID_AutoSaveRaw = wx.NewId()
    
    #Toolbar:
    ID_SaveImage = wx.NewId()
    #ID_Autosave = wx.NewId()
    
    ID_ReloadButton = wx.NewId()
    ID_Autoreload = wx.NewId()
    
    ID_RecordData = wx.NewId()
    
    ID_ImagingChoice = wx.NewId()
    ID_Imaging_AB = wx.NewId()
    ID_ExpansionTimeK = wx.NewId()
    ID_ExpansionTimeNa = wx.NewId()
    ID_OpticalDensityMaxK = wx.NewId()
    ID_OpticalDensityMaxNa = wx.NewId()
    
    ID_DifferenceK= wx.NewId()
    ID_DifferenceNa = wx.NewId()
    ID_FilterK = wx.NewId()
    ID_FilterNa = wx.NewId()
    ID_FilterScalarK = wx.NewId()
    ID_FilterScalarNa = wx.NewId()
    ID_FilterSigmaK = wx.NewId()
    ID_FilterSigmaNa = wx.NewId()
    ID_SaveBitmap = wx.NewId()
    ID_ReferenceK = wx.NewId()
    ID_ReferenceNa = wx.NewId()
    ID_CreateReference = wx.NewId()
    
    ID_CMapChoiceK = wx.NewId()
    ID_CMapChoiceNa = wx.NewId()
    
    ID_RotateK = wx.NewId()
    ID_RotateNa = wx.NewId()
    
    ID_VminK = wx.NewId()
    ID_VminNa = wx.NewId()
    
    ID_VmaxK = wx.NewId()
    ID_VmaxNa = wx.NewId()
    
    ID_AutoscaleClim = wx.NewId()
    
    #ImageTree
    ID_AnalyzeImageButton = wx.NewId()
    ID_ImageTreeReloadSavedImage = wx.NewId()
    ID_ImageTreeRescan = wx.NewId()
    #
    
    #other settings
    imagefilename = settings.imagefile
    rawimg1filename = settings.rawimage1file
    rawimg2filename = settings.rawimage2file
    rawimg3filename = settings.rawimage3file
    autosave_abs = False
    autosave_raw = False

    
    roiK = [ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900)] #list of ROIs for all imaging settings
    roiNa = [ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900)]



    imaging_parlist = [{'K': imagingpars.ImagingParsHorizontalHRNa(), 'Na': imagingpars.ImagingParsHorizontalNa()},
                       {'K': imagingpars.ImagingParsHorizontalNa(), 'Na': imagingpars.ImagingParsHorizontalHRNa()},

                       {'K': imagingpars.ImagingParsAxialNa(), 'Na': imagingpars.ImagingParsHorizontalNa()},
                       {'K': imagingpars.ImagingParsAxialNa(), 'Na': imagingpars.ImagingParsHorizontalHRNa()},
                       {'Na': imagingpars.ImagingParsAxialNa(), 'K': imagingpars.ImagingParsHorizontalNa()},
                       {'Na': imagingpars.ImagingParsAxialNa(), 'K': imagingpars.ImagingParsHorizontalHRNa()},

                       {'K': imagingpars.ImagingParsCMOS(), 'Na': imagingpars.ImagingParsHorizontalNa()},
                       {'K': imagingpars.ImagingParsHorizontalNa(), 'Na': imagingpars.ImagingParsCMOS()},
                       
                       {'K': imagingpars.ImagingParsVerticalNa(), 'Na': imagingpars.ImagingParsHorizontalNa()},
                       {'Na': imagingpars.ImagingParsVerticalNa(), 'K': imagingpars.ImagingParsHorizontalNa()},

                       {'K': imagingpars.ImagingParsAxialNa(), 'Na': imagingpars.ImagingParsVerticalNa()},
                       {'K': imagingpars.ImagingParsVerticalNa(), 'Na': imagingpars.ImagingParsAxialNa()},

                       {'Na': imagingpars.ImagingParsHorizontalNa(), 'K': imagingpars.ImagingParsHorizontalNa()},
                       {'Na': imagingpars.ImagingParsVerticalNa(), 'K': imagingpars.ImagingParsVerticalNa()},

#                       {'K': imagingpars.ImagingParsVerticalK(), 'Na': imagingpars.ImagingParsHorizontalK()},
#                       {'Na': imagingpars.ImagingParsVerticalK(), 'K': imagingpars.ImagingParsHorizontalK()},
#                       {'K': imagingpars.ImagingParsHorizontalK(), 'Na': imagingpars.ImagingParsHorizontalK()},
#                       {'K': imagingpars.ImagingParsVerticalK(), 'Na': imagingpars.ImagingParsVerticalK()},
#                       
#                       {'Na': imagingpars.ImagingParsVerticalNa(), 'K': imagingpars.ImagingParsHorizontalK()},
#                       {'Na': imagingpars.ImagingParsVerticalK(), 'K': imagingpars.ImagingParsHorizontalNa()},
#                       {'K': imagingpars.ImagingParsVerticalNa(), 'Na': imagingpars.ImagingParsHorizontalK()},
#                       {'K': imagingpars.ImagingParsVerticalK(), 'Na': imagingpars.ImagingParsHorizontalNa()},
                        ]


#-------------------------
  
    imaging_pars = imaging_parlist[0]
    markersNapositions = []
    markersKpositions = []
    


    def select_imaging(self, n):
        """
        select new imaging parameters. sync controls with new settings. 
        keep imaging pars of fitting objects in sync
        """
        
        self.imaging_pars = self.imaging_parlist[n]
        self.sync_imaging_pars_expansion_time(updateControls=True)
        self.sync_imaging_pars_maximum_optical_density(updateControls=True)

        self.Na.fit.set_imaging_pars(self.imaging_pars['Na'])
        self.Na.fitpars.imaging_pars = self.imaging_pars['Na']
        
        self.K.fit.set_imaging_pars(self.imaging_pars['K'])
        self.K.fitpars.imaging_pars = self.imaging_pars['K']

    def sync_imaging_pars_expansion_time(self, updateControls=False):
        "keep imaging_pars in sync with controls for expansion time"
        if updateControls:
            self.expansion_time_K.Value = self.imaging_pars['K'].expansion_time
            self.expansion_time_Na.Value = self.imaging_pars['Na'].expansion_time
        else:
            self.imaging_pars['K'].expansion_time = self.expansion_time_K.Value
            self.imaging_pars['Na'].expansion_time = self.expansion_time_Na.Value
        
    def sync_imaging_pars_maximum_optical_density(self, updateControls=False):
        if updateControls == False:
            #get values from control
            try:
                ODmaxK = float(self.entry_optical_density_max_K.Value)
            except ValueError:
                ODmaxK = 0
                
            try:
                ODmaxNa = float(self.entry_optical_density_max_Na.Value)
            except ValueError:
                ODmaxNa = 0
            
            self.imaging_pars['K'].ODmax = ODmaxK
            self.imaging_pars['Na'].ODmax = ODmaxNa
        else:
            #update controls
            ODmaxK = self.imaging_pars['K'].ODmax
            ODmaxNa = self.imaging_pars['Na'].ODmax
        
        #always update control
        self.entry_optical_density_max_K.Value = "%.1f" % ODmaxK
        self.entry_optical_density_max_Na.Value = "%.1f" % ODmaxNa    

    def select_roi(self, sel):
        self.Na.set_roi(self.roiNa[sel])
        self.K.set_roi(self.roiK[sel])
    
    def OnChoiceCMap(self,event):
        if event.Id==self.ID_CMapChoiceK:
            sel=self.cmap_choiceK.GetStringSelection()
            self.imaging_pars['K'].palette=sel
            self.K.himg.set_cmap(pylab.get_cmap(self.imaging_pars['K'].palette))
            self.K.update_image()
            self.K.redraw()
            
        elif event.Id==self.ID_CMapChoiceNa:
            sel=self.cmap_choiceNa.GetStringSelection()
            self.imaging_pars['Na'].palette=sel
            self.Na.himg.set_cmap(pylab.get_cmap(self.imaging_pars['Na'].palette))
            self.Na.update_image()
            self.Na.redraw()
        
    def OnChoiceImaging(self, event):
        sel = self.imaging_choice.GetSelection()

        self.select_imaging(sel)
        self.select_roi(sel)

        self.Na.center_roi()
        self.K.center_roi()

        self.Na.update()
        self.K.update()
        
        self.Na.himg.set_cmap(pylab.get_cmap(self.imaging_pars['Na'].palette))
        self.K.himg.set_cmap(pylab.get_cmap(self.imaging_pars['K'].palette))
        self.cmap_choiceK.SetStringSelection(self.imaging_pars['K'].palette)
        self.cmap_choiceNa.SetStringSelection(self.imaging_pars['Na'].palette)

        
    def OnChoice_AB(self, event):
        self.load_image(self.imagefilename)
        
    def OnChangeExpansionTime(self, event):
        id = event.GetId()
        self.sync_imaging_pars_expansion_time()
        if id == self.ID_ExpansionTimeK:
            self.K.update_parameters()
            self.K.redraw()
            self.results.UpdateResults({'K': self.K.fitpars})
        elif id == self.ID_ExpansionTimeNa:
            self.Na.update_parameters()
            self.Na.redraw()
            self.results.UpdateResults({'Na': self.Na.fitpars})
        
    def OnChangeMaximumOpticalDensity(self, event):
        id = event.Id
        self.sync_imaging_pars_maximum_optical_density()
        #TODO: need reload image
        if id == self.ID_OpticalDensityMaxK:
            self.K.update()
        elif id == self.ID_OpticalDensityMaxNa:
            self.Na.update()


    def OnChangeRotation(self, event):
        id = event.Id
        #TODO: need reload image
        if id == self.ID_RotateK:
            self.K.rotation_angle = float(self.rotation_K.GetValue())
            if self.K.rotation_angle != 0.0:
                self.K.rotation_active = True
            else:
                self.K.rotation_active = False
            self.K.update()
        elif id == self.ID_RotateNa:
            self.Na.rotation_angle = float(self.rotation_Na.GetValue())
            if self.Na.rotation_angle != 0.0:
                self.Na.rotation_active = True
            else:
                self.Na.rotation_active = False
            self.Na.update()
            
            
    def OnChangeVmin(self, event):
        id = event.Id
        #TODO: need reload image
        if id == self.ID_VminK:
            vmin = float(self.vmin_K.GetValue())
            self.K.himg.set_clim(vmin=vmin)
            self.K.redraw()
        elif id == self.ID_VminNa:
            vmin = float(self.vmin_Na.GetValue())
            self.Na.himg.set_clim(vmin=vmin)
            self.Na.redraw()
            
    def OnChangeVmax(self, event):
        id = event.Id
        #TODO: need reload image
        if id == self.ID_VmaxK:
            vmax = float(self.vmax_K.GetValue())
            self.K.himg.set_clim(vmax=vmax)
            self.K.redraw()
        elif id == self.ID_VmaxNa:
            vmax = float(self.vmax_Na.GetValue())
            self.Na.himg.set_clim(vmax=vmax)
            self.Na.redraw()
            
    def OnAutoscaleClim(self, event):
        id = event.Id
        print 'auto clim'
        if id == self.ID_AutoscaleClim:
            for window in [self.K, self.Na]:
                im = window.himg.get_array()
                m, s = im.mean(), im.std()
                window.himg.set_clim(m-s, m+s)
                print window.himg.get_clim()
                window.redraw()


    def OnChangeFilterSigma(self, event):
        id = event.Id
        if id == self.ID_FilterSigmaK:
            self.K.filter_sigma=float(self.entry_filter_sigma_K.GetValue())
            self.K.update_image()
            self.K.redraw()
        elif id == self.ID_FilterSigmaNa:
            self.Na.filter_sigma=float(self.entry_filter_sigma_Na.GetValue())
            self.Na.update_image()
            self.Na.redraw()
        self.UpdateResiduals()
    
    def UpdateResiduals(self):
        try:
            if self.K.target_results!=None:
                row=self.K.target_results["row"]
            else:
                row=self.results.active_row
            
            if(self.Na.bec_residuals!= None):    
                self.results.SetValueNamed(row, "resid Na", self.Na.bec_residuals)
            
            if(self.K.bec_residuals!= None):    
                self.results.SetValueNamed(row, "resid K", self.K.bec_residuals)
            
            self.results.GetView().Refresh()
        except AttributeError:
            print "error while updating residuals"
            
    def OnChangeFilterScalar(self, event):
        id = event.Id
        if id == self.ID_FilterScalarK:
            self.K.filter_scalar=float(self.entry_filter_scalar_K.GetValue())
            self.K.update_image()
            self.K.redraw()
        elif id == self.ID_FilterScalarNa:
            self.Na.filter_scalar=float(self.entry_filter_scalar_Na.GetValue())
            self.Na.update_image()
            self.Na.redraw()
        #self.UpdateResiduals()
            
    def OnMenuResultsSave(self, event):
        id = event.GetId()

        filename = self.results.filename
        if id == self.ID_ResultsSaveAs or filename is None:
            filename = self.results_save_filename_ask(self.results)
        
        self.results_save(filename, self.results)
        
    def results_autosave(self):
        filename = self.results.filename
        print "autosaving %s results"%filename    
        if filename is None:
            print "skipping"
        else:
            self.results_save_overwrite(filename, self.results)

    def results_save_filename_ask(self, results):
        """return filename or None"""
        # search for date in results.name
        expr = re.compile(r'(?P<day>[0-9]{8}-|[0-9]{4}-[0-9]{2}-[0-9]{2}-)?(?P<name>.*)')
        print "name: %s"%self.results.name
        m = expr.match(self.results.name)
        print '<Match: %r, groups=%r>' % (m.group(), m.groups())
        if m is not None and m.group('day') is not None:    
            day = ''
        else:
            day = time.strftime("%Y-%m-%d") + '-'
        filename_proposal =  day + self.results.name 
        savedialog = wx.FileDialog(
            self.frame,
            message="Save results as ...",
            defaultDir=self.get_data_dir(),
            defaultFile=filename_proposal,
            wildcard="CSV file (*.csv)|*.csv|All files (*.*)|*.*",
            style=wx.SAVE | wx.FD_CHANGE_DIR)

        if savedialog.ShowModal() == wx.ID_OK:
            filename = savedialog.GetPath()
        else:
            filename = None

        savedialog.Destroy()
        return filename

    def results_save(self, filename, results, forcesave = False):#False):

        """save results to file, check for overwriting. If forcesave
        is True, don't give up."""

        success = False
        while not success:

            if filename is None:
                if forcesave:
                    filename = self.results_save_filename_ask(results)
                    if filename is None:
                        #continue
                        print "data not saved"
                        break
                else:
                    print "data not saved!"
                    results.filename = None
                    return success

            #check if file exists
            if filename is not None and os.access(filename, os.F_OK):
                #file already exists, ask to overwrite it
                MB = wx.MessageDialog(self.frame,
                                      "File " + filename + 
                                      " already exists. \nDo you want to overwrite it?",
                                      caption="Save File ...",
                                      style=wx.YES_NO | wx.ICON_EXCLAMATION,
                                      )
                answer = MB.ShowModal()
                MB.Destroy()
                if answer == wx.ID_YES:
                    pass
                elif forcesave:
                    #don't overwrite, ask again
                    filename = None
                    continue
                else:
                    print "data file not saved"
                    results.filename = None
                    return success
                
            try:
                results.save_data_csv(filename)
            except StandardError, e:
                import traceback
                traceback.print_exc()
                msg = wx.MessageDialog(self.frame,
                                       'Error saving data to file: "%s"\n%s'\
                                       %(filename, traceback.format_exc(e)),
                                       style = wx.OK | wx.ICON_ERROR,
                                       )
                msg.ShowModal()
                success = False
                filename = None
                results.filename = None
            else:
                results.filename = filename
                success = True

            return success
            
    def results_save_overwrite(self, filename, results,):
        """save results to file, overwriting."""
        success = False
        while not success:
            #check if file exists
            if filename is not None and os.access(filename, os.F_OK):
                #file already exists,
                print "overwriting %s"%filename
            try:
                results.save_data_csv(filename)
            except StandardError, e:
                import traceback
                traceback.print_exc()
                msg = wx.MessageDialog(self.frame,
                                       'Error saving data to file: "%s"\n%s'\
                                       %(filename, traceback.format_exc(e)),
                                       style = wx.OK | wx.ICON_ERROR,
                                       )
                msg.ShowModal()
                success = False
                filename = None
                results.filename = None
            else:
                results.filename = filename
                success = True
                
            return success
                
    def OnMenuResultsLoad(self, event):
        loaddialog = wx.FileDialog(
            self.frame,
            message="Load results",
            defaultDir=settings.imagesavepath,
            defaultFile='',
            wildcard="CSV file (*.csv)|*.csv",
            style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST | wx.FD_CHANGE_DIR)
        if loaddialog.ShowModal() == wx.ID_OK:
            fullpath = loaddialog.GetPath()
            dir, filename = os.path.split(fullpath)
            name, ext = os.path.splitext(filename)
            self.select_measurement(name)
            self.results.load_data_csv(fullpath)

        loaddialog.Destroy()

    def OnMenuResultsSaveTemplate(self, event):
        savedialog = wx.FileDialog(
            self.frame,
            message="Save settings to template file ...",
            defaultDir=settings.templatedir,
            defaultFile='template' + self.results.name,
            wildcard="template file (*.tpl)|*.tpl|All files (*.*)|*.*",
            style=wx.SAVE)

        if savedialog.ShowModal() == wx.ID_OK:
            templatefile = savedialog.GetPath()

            #save metadata to template file
            metadata = self.results.give_metadata()
            output = open(templatefile, 'wb')
            pickle.dump(metadata, output)
            output.close()

        savedialog.Destroy()


    def OnMenuResultsApplyTemplate(self, event):
        loaddialog = wx.FileDialog(
            self.frame,
            message="Load settings from template file ...",
            defaultDir=settings.templatedir,
            defaultFile='',
            wildcard="template file (*.tpl)|*.tpl|All files (*.*)|*.*",
            style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)

        if loaddialog.ShowModal() == wx.ID_OK:
            infile = open(loaddialog.GetPath(), 'rb')
            metadata = pickle.load(infile)
            infile.close()
            self.results.GetView().ApplyTemplate(metadata)

        loaddialog.Destroy()

    def OnMenuResultsNew(self, event):
        """
        Create new measurement, based on actual measurement (copy settings).
        """

        #get metadata from active measurement
        metadata = self.results.give_metadata()

        #TODO: check if saved, save

        #make proposal for measurement
        allnames = self.measurement_combobox.Strings
        name = self.results.name #name of active measurement

        #try to match name like "meas 5"
        pattern = r"(?P<name>[a-zA-Z_]+\s*)(?P<nr>\d+)"
        match = re.match(pattern, name)
        if match is not None:
            #if name matches pattern "meas 5", use "meas 6" as
            #proposal. If this already used, try higher numbers
            gd = match.groupdict()
            prefix = gd['name']
            nr = int(gd['nr'])
            for k in range(100):
                proposal = "%s%d" % (prefix, nr + k)
                if proposal not in allnames:
                    break
        else:
            for k in range(2, 100):
                proposal = "%s %d" % (name.strip(), k)
                if proposal not in allnames:
                    break
        
        #ask for new measurement name, with proposal
        dlg = wx.TextEntryDialog(
            self.frame,
            message="New measurement name:",
            caption="Create new measurement ...",
            defaultValue=proposal)
        if dlg.ShowModal() == wx.ID_OK:
            newname = dlg.Value
        else:
            newname = None
        dlg.Destroy()

        #create new measurement (update combo box)
        #plot results

        if newname:
            self.select_measurement(newname)
            self.results.GetView().ApplyTemplate(metadata)
            self.plot_results()
    
    def OnCreateReferenceFromMask(self,event):
        table=self.gridpanel.pages[self.gridpanel.notebook.GetSelection()].Table
        maskedrows=table.rowmask
        maskedrowsfiles= [table.GetValueNamed(i, 'Filename') for i in range(len(maskedrows)) if maskedrows[i]]
        imagesavedir = self.get_data_dir(subdir = 'images')
        maskedrowsfilespath=[]
        for row in range(len(maskedrowsfiles)):
            if(maskedrowsfiles[row]!=""):
                maskedrowsfilespath.append(os.path.normpath(os.path.join(imagesavedir, maskedrowsfiles[row])))
        
        if len(maskedrowsfilespath)>0:
            mean_sis.create_reference(maskedrowsfilespath)
        
        
    def OnCreateReference(self, event):
        loaddialog = wx.FileDialog(
            self.frame,
            message="Select references",
            defaultDir=settings.imagesavepath,
            defaultFile='',
            wildcard="SIS file (*.sis)|*.sis",
            style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST | wx.FD_CHANGE_DIR | wx.FD_MULTIPLE)
        if loaddialog.ShowModal() == wx.ID_OK:
            fullpath = loaddialog.GetPaths()
            mean_sis.create_reference(fullpath)

        loaddialog.Destroy()    
        
    def OnReferenceSelect(self,event):
        if event.Id==self.ID_ReferenceNa:
            if event.IsChecked():
                self.Na.reference_diff=True
                _, self.Na.reference_img = loadimg(settings.referencefile)
            else:
                self.Na.reference_diff=False
            
            self.Na.update_image()
            self.Na.redraw()
            
        elif event.Id==self.ID_ReferenceK:
            if event.IsChecked():
                self.K.reference_diff=True
                self.K.reference_img, _= loadimg(settings.referencefile)
            else:
                self.K.reference_diff=False
            
            self.K.update_image()
            self.K.redraw()

                   
    def OnDifferenceSelectNa(self,event):
        if event.IsChecked():
            self.Na.difference=True
        else:
            self.Na.difference=False
        self.Na.update()
            
    def OnDifferenceSelectK(self,event):
        if event.IsChecked():
            self.K.difference=True
        else:
            self.K.difference=False
        self.K.update()
        
    def OnFilterSelectK(self,event):
        if event.IsChecked():
            self.K.filter=True
        else:
            self.K.filter=False
        self.K.update_image()
        self.K.redraw()
        
        self.UpdateResiduals()
        
    def OnFilterSelectNa(self,event):
        if event.IsChecked():
            self.Na.filter=True
        else:
            self.Na.filter=False
        self.Na.update_image()
        self.Na.redraw()
        
        self.UpdateResiduals()

    def OnMenuShow(self, event):
        id = event.GetId()

        if id in [self.ID_ShowNa]:
            pane = self.mgr.GetPane("Na")
            if not pane.IsShown():
                pane.Float().Show()
                self.mgr.Update()

        if id in [self.ID_ShowK]:
            pane = self.mgr.GetPane("K")
            if not pane.IsShown():
                pane.Float().Show()
                self.mgr.Update()

        if id in [self.ID_ShowTools]:
            pane = self.mgr.GetPane('toolbar')
            if not pane.IsShown():
                pane.Show().ToolbarPane()
                self.mgr.Update()

        if id in [self.ID_ShowSavedImages]:
            panetree = self.mgr.GetPane('imagetree')
            paneimgna = self.mgr.GetPane('savedimageNa')
            paneimgk = self.mgr.GetPane('savedimageK')
            if not panetree.IsShown():
                panetree.Right().Show()

            if not paneimgna.IsShown():
                paneimgna.Float().Show()

            if not paneimgk.IsShown():
                paneimgk.Float().Show()

            self.mgr.Update()
            
    def OnMenuPerspective(self, event):
        #print self.mgr.SavePerspective()
        id = event.GetId()
        if id in [self.ID_PerspectiveFull]:
            self.mgr.LoadPerspective(self.perspective_full)
        if id in [self.ID_PerspectiveNa]:
            self.mgr.LoadPerspective(self.perspective_Na)
        if id in [self.ID_PerspectiveK]:
            self.mgr.LoadPerspective(self.perspective_K)

    def OnMenuFit(self, event):
        id = event.GetId()
        self.markersNapositions = [self.Na.markers[0].position[0],
                          self.Na.markers[0].position[1]]
        self.markersKpositions = [self.K.markers[0].position[0],
                          self.K.markers[0].position[1]]
        # side markers Na
        for k, marker in enumerate(self.Na.markers):
            if k > 0:
                self.markersNapositions.append(self.Na.markers[k].position[0])
                self.markersNapositions.append(self.Na.markers[k].position[1])
                # insert automatically symmetric side marker
                # TO DO: check if out of range
                self.markersNapositions.append(2*self.markersNapositions[0]-self.Na.markers[k].position[0])
                self.markersNapositions.append(2*self.markersNapositions[1]-self.Na.markers[k].position[1])
        # side markers K
        for k, marker in enumerate(self.K.markers):
            if k > 0: 
                self.markersKpositions.append(self.K.markers[k].position[0])
                self.markersKpositions.append(self.K.markers[k].position[1])
                # insert automatically symmetric side marker
                # TO DO: check if out of range
                self.markersKpositions.append(2*self.markersKpositions[0]-self.K.markers[k].position[0])
                self.markersKpositions.append(2*self.markersKpositions[1]-self.K.markers[k].position[1])
        
                
        print 'markers Na', self.markersNapositions
        print 'markers K', self.markersKpositions
            

        if id in [self.ID_FitShowContours]:
            self.Na.show_contours = event.IsChecked()
            self.K.show_contours = event.IsChecked()
            self.Na.update()
            self.K.update()

        if id in [self.ID_FitNaNone]:
            self.Na.fit = fitting.NoFit(self.imaging_pars['Na'])
            self.Na.update()
            
        if id in [self.ID_NaImageIntegration]:
            self.Na.fit = fitting.ImageIntegration(self.imaging_pars['Na'])
            self.Na.update()
            
        if id in [self.ID_FitNaGauss1D]:
            self.Na.fit = fitting.Gauss1d(self.imaging_pars['Na'])
            self.Na.update()

        if id in [self.ID_FitNaLorentzGauss1D]:
            self.Na.fit = fitting.LorentzGauss1d(self.imaging_pars['Na'])
            self.Na.update()

        if id in [self.ID_FitNaGaussGauss1D]:
            self.Na.fit = fitting.GaussGauss1d(
                        self.imaging_pars['Na'], self.markersNapositions)
            self.Na.update()

        if id in [self.ID_FitNaGauss]:
            self.Na.fit = fitting.Gauss2d(self.imaging_pars['Na'])
            self.Na.update()

        if id in [self.ID_FitNaGaussSym]:
            self.Na.fit = fitting.GaussSym2d(self.imaging_pars['Na'])
            self.Na.update()

        if id in [self.ID_FitNaGaussBose]:
            self.Na.fit = fitting.GaussBose2d(self.imaging_pars['Na'])
            self.Na.update()

        if id in [self.ID_FitNaTF]:
            self.Na.fit = fitting.ThomasFermi2d(self.imaging_pars['Na'])
            self.Na.update()     
            
        if id in [self.ID_FitNaTFSlice]:
            self.Na.fit = fitting.ThomasFermiSlice(self.imaging_pars['Na'])
            self.Na.update()     

        if id in [self.ID_FitNaBimodal]:
            self.Na.fit = fitting.Bimodal2d(self.imaging_pars['Na'])
            self.Na.update()
        
        if id in [self.ID_FitNaBimodalSplit]:
            self.Na.fit = fitting.Bimodal2dSplit(self.imaging_pars['Na'])
            self.Na.update()
            
        if id in [self.ID_FitNaThomasFermiHole]:
            self.Na.fit = fitting.ThomasFermiHole2d(self.imaging_pars['Na'])
            self.Na.update()
            
        if id in [self.ID_FitNaBimodalIntegrate1dSplit]:
            self.Na.fit = fitting.BimodalIntegrate1dSplit(self.imaging_pars['Na'])
            self.Na.update()
        
        if id in [self.ID_FitNaGaussIntegrate1dPinning]:
            self.Na.fit = fitting.GaussIntegrate1dPinning(self.imaging_pars['Na'])
            self.Na.update()

        if id in [self.ID_FitNaBimodalGaussGauss]:
            state = self.menu.FindItemById(self.ID_MarkerGuessNa).IsChecked()
            self.Na.set_marker_guess(state)
            if self.Na.marker_guess:
                self.Na.fit = fitting.BimodalGaussGauss2d(
                            self.imaging_pars['Na'], self.markersNapositions
                            )
            else:
                self.Na.fit = fitting.BimodalGaussGauss2d(
                    self.imaging_pars['Na'])
            self.Na.update()

        if id in [self.ID_FitNaBoseBimodal]:
            self.Na.fit = fitting.BoseBimodal2d(self.imaging_pars['Na'])
            self.Na.update()

        if id in [self.ID_FitNaNGauss]:
            self.Na.fit = fitting.BimodalGaussGauss2d(
                            self.imaging_pars['Na'], self.markersNapositions
                            )
            self.Na.fit = fitting.NGauss2d(self.imaging_pars['Na'], 
                            self.markersNapositions
                            )
            self.Na.update()

        if id in [self.ID_FitKNone]:
            self.K.fit = fitting.NoFit(self.imaging_pars['K'])
            self.K.update()
            
        if id in [self.ID_KImageIntegration]:
            self.K.fit = fitting.ImageIntegration(self.imaging_pars['K'])
            self.K.update()
            """
            try:
                f = open("test",'w')
                f.write(self.K.fit[2])
                f.close()
            except:
                exit(0)
                """
            
        if id in [self.ID_FitKGauss1D]:
            self.K.fit = fitting.Gauss1d(self.imaging_pars['K'])
            self.K.update()

        if id in [self.ID_FitKLorentzGauss1D]:
            self.K.fit = fitting.LorentzGauss1d(self.imaging_pars['K'], )
            self.K.update()

        if id in [self.ID_FitKGaussGauss1D]:
            self.K.fit = fitting.GaussGauss1d(
                       self.imaging_pars['K'], self.markersKpositions)
            self.K.update()

        if id in [self.ID_FitKGauss]:
            self.K.fit = fitting.Gauss2d(self.imaging_pars['K'])
            self.K.update()

        if id in [self.ID_FitKGaussSym]:
            self.K.fit = fitting.GaussSym2d(self.imaging_pars['K'])
            self.K.update()
            
        if id in [self.ID_FitKGaussBose]:
            self.K.fit = fitting.GaussBose2d(self.imaging_pars['K'])
            self.K.update()

        if id in [self.ID_FitKTF]:
            self.K.fit = fitting.ThomasFermi2d(self.imaging_pars['K'])
            self.K.update()       
        
        if id in [self.ID_FitKTFSlice]:
            self.K.fit = fitting.ThomasFermiSlice(self.imaging_pars['K'])
            self.K.update()   

        if id in [self.ID_FitKBimodal]:
            print 'bimodal K'
            self.K.fit = fitting.Bimodal2d(self.imaging_pars['K'])
            self.K.update()
            
        if id in [self.ID_FitKBimodalSplit]:
            print 'bimodal K - split mxy'
            self.K.fit = fitting.Bimodal2dSplit(self.imaging_pars['K'])
            self.K.update()
            
        if id in [self.ID_FitKThomasFermiHole]:
            print 'bimodal K - split mxy'
            self.K.fit = fitting.ThomasFermiHole2d(self.imaging_pars['K'])
            self.K.update()
            
        if id in [self.ID_FitKBimodalIntegrate1dSplit]:
            self.K.fit = fitting.BimodalIntegrate1dSplit(self.imaging_pars['K'])
            self.K.update()
            
        if id in [self.ID_FitKGaussIntegrate1dPinning]:
            self.K.fit = fitting.GaussIntegrate1dPinning(self.imaging_pars['K'])
            self.K.update()

        if id in [self.ID_FitKBimodalGaussGauss]:
            state = self.menu.FindItemById(self.ID_MarkerGuessK).IsChecked()
            self.K.set_marker_guess(state)
            if self.K.marker_guess:
                self.K.fit = fitting.BimodalGaussGauss2d(
                            self.imaging_pars['K'], self.markersKpositions
                            )
            else:
                self.K.fit = fitting.BimodalGaussGauss2d(
                    self.imaging_pars['K'])
            self.K.update()
        
        if id in [self.ID_FitKBoseBimodal]:
            self.K.fit = fitting.BoseBimodal2d(self.imaging_pars['K'])
            self.K.update()
        


    def OnMenuCompensate(self, event):
        id = event.GetId()

        if id in [self.ID_CompensateOffsetNa]:
            state = self.menu.FindItemById(id).IsChecked()
            self.Na.set_compensate_offset(state)

        if id in [self.ID_CompensateOffsetK]:
            if self.menu.FindItemById(id).IsChecked():
                self.K.set_compensate_offset(True)
            else:
                self.K.set_compensate_offset(False)

    def OnMenuMarkerGuess(self, event):
        id = event.GetId()

        if id in [self.ID_MarkerGuessNa]:
            state = self.menu.FindItemById(id).IsChecked()
            self.Na.set_marker_guess(state)
        if id in [self.ID_MarkerGuessK]:
            if self.menu.FindItemById(id).IsChecked():
                self.K.set_marker_guess(True)
            else:
                self.K.set_marker_guess(False)

    def OnMenuContrast(self, event):
        id = event.GetId()

        K_contrast_dict = {
            self.ID_K_ContrastUltra:  0.2,
            self.ID_K_ContrastHigh:   0.5,
            self.ID_K_ContrastNormal: 1.5,
            self.ID_K_ContrastLow:    2.5,
            }

        Na_contrast_dict = {
            self.ID_Na_ContrastUltra:  0.2,
            self.ID_Na_ContrastHigh:   0.5,
            self.ID_Na_ContrastNormal: 1.5,
            self.ID_Na_ContrastLow:    2.5,
            }

        cmax = K_contrast_dict.get(id)
        if cmax:
            self.K.himg.set_clim(vmax=cmax)
            self.K.redraw()
            #print id, cmax
            
        cmax = Na_contrast_dict.get(id)
        if cmax:
            self.Na.himg.set_clim(vmax=cmax)
            self.Na.redraw()

    def OnMenuMarker(self, event):
        markerindex = [self.ID_MarkerA, self.ID_MarkerB, self.ID_MarkerC].index(event.Id)
        self.Na.markers[markerindex].set_visible(event.IsChecked())
        self.K.markers[markerindex].set_visible(event.IsChecked())
        self.Na.redraw()
        self.K.redraw()

    def OnMenuReload(self, event):
        id = event.GetId()

        if id in [self.ID_Reload]:
            self.PostReloadImageEvent()

    def OnMenuAutoSave(self, event):
        id = event.GetId()

        if id in [self.ID_AutoSaveAbs]:
            self.autosave_abs = event.IsChecked()

        if id in [self.ID_AutoSaveRaw]:
            self.autosave_raw = event.IsChecked()

    def OnReloadButton(self, event):
        self.PostReloadImageEvent()

    def OnSaveImage(self, event):
        evt = SaveImageEvent()
        wx.PostEvent(self.frame, evt)

###########
    def OnAutosaveToggle(self, event):
        """Callback for Auto Save Image Checkbox. Stores status to
        self.autosave"""
        if event.IsChecked():
            self.autosave = True
        else:
            self.autosave = False
###############

    def get_data_dir(self, subdir = ''):
        """change current directory to data dir"""
        directory = os.path.join(settings.imagesavepath,
                                 time.strftime("%Y/%Y-%m-%d/"),
                                 subdir)
        if not os.access(directory, os.F_OK):
            try:#try to create dir
                os.makedirs(directory)
            except OSError:
                print "cannot create data dir"
                return os.getcwd()
        return directory
    
    def OnSaveBitmap(self, event):
        imagesavedir = self.get_data_dir(subdir = 'bitmaps')
        if self.K.target_results!=None:
            name=self.K.target_results["name"]
            row=self.K.target_results["row"]
        else:
            row=self.results.active_row
            name=self.results.name
        
        imagesavefilename0 = "%s-%s-%04d_0.png" % (time.strftime("%Y-%m-%d-T%H%M%S"),
                                     name,
                                     int(row))
        imagesavefilename1 = "%s-%s-%04d_1.png" % (time.strftime("%Y-%m-%d-T%H%M%S"),
                                     name,
                                     int(row))
        
        
        imagesavefilename01 = "%s-%s-%04d.sis" % (time.strftime("%Y-%m-%d-T%H%M%S"),
                                     name,
                                     int(row))
        
        imagesavefilenamefull01 = os.path.normpath(os.path.join(imagesavedir, imagesavefilename01))
        
        
        imagesavefilenamefull0 = os.path.normpath(os.path.join(imagesavedir, imagesavefilename0))
        imagesavefilenamefull1 = os.path.normpath(os.path.join(imagesavedir, imagesavefilename1))
        
        #test if file already exists
        if os.access(imagesavefilenamefull0, os.F_OK):
            MB = wx.MessageDialog(self.frame,
                                  "Image file " + imagesavefilename0 + 
                                  " already exists. \nDo you want to overwrite it?",
                                  caption="Save Image File ...",
                                  style=wx.YES_NO | wx.ICON_EXCLAMATION,

                                  )
            answer = MB.ShowModal()
            MB.Destroy()
            if  answer == wx.ID_YES:
                pass
            else:
                print "image file not saved"
                return
        if os.access(imagesavefilenamefull1, os.F_OK):
            MB = wx.MessageDialog(self.frame,
                                  "Image file " + imagesavefilename1 + 
                                  " already exists. \nDo you want to overwrite it?",
                                  caption="Save Image File ...",
                                  style=wx.YES_NO | wx.ICON_EXCLAMATION,

                                  )
            answer = MB.ShowModal()
            MB.Destroy()
            if  answer == wx.ID_YES:
                pass
            else:
                print "image file not saved"
                return

        
        if os.access(imagesavefilenamefull01, os.F_OK):
            MB = wx.MessageDialog(self.frame,
                                  "Image file " + imagesavefilename01 + 
                                  " already exists. \nDo you want to overwrite it?",
                                  caption="Save Image File ...",
                                  style=wx.YES_NO | wx.ICON_EXCLAMATION,

                                  )
            answer = MB.ShowModal()
            MB.Destroy()
            if  answer == wx.ID_YES:
                pass
            else:
                print "image file not saved"
                return
        

        K_bit=self.K.create_bitmap()
        Na_bit=self.Na.create_bitmap()

        NaK_bit=numpy.concatenate((numpy.ma.filled(K_bit,0), numpy.ma.filled(Na_bit,0)))
        NaK_bit=(NaK_bit+1.0)*6553.6

        write_raw_image(imagesavefilenamefull01,NaK_bit)
        
        self.K.save_bitmap(imagesavefilenamefull0,pylab.get_cmap(self.imaging_pars['K'].palette))
        self.Na.save_bitmap(imagesavefilenamefull1,pylab.get_cmap(self.imaging_pars['Na'].palette))
            
    def OnSaveImageEvent(self, event):
        imagesavedir = self.get_data_dir(subdir = 'images')
        ts_full = time.strftime("%Y-%m-%d-T%H%M%S")
        ts_time = time.strftime("%H:%M:%S")
        imagesavefilename = "%s-%s-%04d.sis" % (ts_full,
                                             self.results.name,
                                             self.results.active_row)
        rawimg1savefilename = "r1-" + imagesavefilename
        rawimg2savefilename = "r2-" + imagesavefilename
        rawimg3savefilename = "r3-" + imagesavefilename

        imagesavefilenamefull = os.path.normpath(os.path.join(imagesavedir, imagesavefilename))
        rawimg1savefilenamefull = os.path.normpath(os.path.join(imagesavedir, rawimg1savefilename))
        rawimg2savefilenamefull = os.path.normpath(os.path.join(imagesavedir, rawimg2savefilename))
        rawimg3savefilenamefull = os.path.normpath(os.path.join(imagesavedir, rawimg3savefilename))

        #test if file already exists
        if os.access(imagesavefilenamefull, os.F_OK):
            MB = wx.MessageDialog(self.frame,
                                  "Image file " + imagesavefilename + 
                                  " already exists. \nDo you want to overwrite it?",
                                  caption="Save Image File ...",
                                  style=wx.YES_NO | wx.ICON_EXCLAMATION,

                                  )
            answer = MB.ShowModal()
            MB.Destroy()
            if  answer == wx.ID_YES:
                pass
            else:
                print "image file not saved"
                return
                                
        shutil.copy2(self.imagefilename, imagesavefilenamefull)
        statusbar_msg = "image saved as: " + imagesavefilename
        if self.autosave_raw:   # save also raw images
            statusbar_msg = statusbar_msg + "  (+ raw images as: r1-..., r2-..., r3-...)"
            shutil.copy2(self.rawimg1filename, rawimg1savefilenamefull)
            shutil.copy2(self.rawimg2filename, rawimg2savefilenamefull) 
            shutil.copy2(self.rawimg3filename, rawimg3savefilenamefull)

        self.savebutton.SetBackgroundColour(wx.NamedColor("GREEN"))
        self.savebutton.Refresh()

        #self.saved_image_name.SetLabel(imagesavefilename)
        self.statusbar.SetStatusText(statusbar_msg)
        self.UpdateResultsFilename(imagesavefilename)
        # added timestamp
        self.UpdateResultsTimestamp(ts_time)

    def OnInit(self):
        ## create main frame

        splash = CamSplashScreen()
        splash.Show()

        self.SetAppName("Cam")
        
        self.frame = wx.Frame(None,
                              title="Image Display and Fit",
                              size=(1300, 750))

        #set main icon
        icons = wx.IconBundle()
        for icon in ['cam16.png',
                     'cam24.png',
                     'cam32.png',
                     'cam48.png',
                     ]:
                   
            icons.AddIconFromFile(os.path.join(settings.bitmappath, icon),
                                  wx.BITMAP_TYPE_PNG)
        self.frame.SetIcons(icons)

        #create manager
        self.mgr = wx.aui.AuiManager(self.frame,
                                     wx.aui.AUI_MGR_RECTANGLE_HINT 
 | wx.aui.AUI_MGR_ALLOW_FLOATING
                                     )

        self.mgr.SetDockSizeConstraint(0.5, 0.75)
        
        ## create menu bar
        
        #submenus Show
        view_menu = wx.Menu()
        view_menu.Append(self.ID_ShowNa, "Show Na")
        view_menu.Append(self.ID_ShowK, "Show K")
        view_menu.Append(self.ID_ShowTools, "Show Toolbar")
        view_menu.Append(self.ID_ShowSavedImages, "Browse saved images")
        self.frame.Bind(wx.EVT_MENU, self.OnMenuShow, id=self.ID_ShowNa)
        self.frame.Bind(wx.EVT_MENU, self.OnMenuShow, id=self.ID_ShowK)
        self.frame.Bind(wx.EVT_MENU, self.OnMenuShow, id=self.ID_ShowTools)
        self.frame.Bind(wx.EVT_MENU, self.OnMenuShow, id=self.ID_ShowSavedImages)

        #subsubmenu Perspectives
        perspective_menu = wx.Menu()
        perspective_menu.Append(self.ID_PerspectiveFull, "Show all")
        perspective_menu.Append(self.ID_PerspectiveNa, "Na only")
        perspective_menu.Append(self.ID_PerspectiveK, "K only")
        self.frame.Bind(wx.EVT_MENU, self.OnMenuPerspective, id=self.ID_PerspectiveFull)
        self.frame.Bind(wx.EVT_MENU, self.OnMenuPerspective, id=self.ID_PerspectiveNa)
        self.frame.Bind(wx.EVT_MENU, self.OnMenuPerspective, id=self.ID_PerspectiveK)

        view_menu.AppendMenu(self.ID_Perspective, "Perspectives", perspective_menu)

        #submenu Display
        display_menu = wx.Menu()
        #subsubmenu K display
        display_K_menu = wx.Menu()
        display_K_menu.AppendRadioItem(self.ID_K_ContrastHigh, 'high', 'high contrast')
        display_K_menu.AppendRadioItem(self.ID_K_ContrastNormal, 'normal', 'low contrast')
        display_K_menu.AppendRadioItem(self.ID_K_ContrastLow, 'low', 'high contrast')
        display_K_menu.AppendRadioItem(self.ID_K_ContrastUltra, 'ultra', 'ultra high contrast')
        display_menu.AppendMenu(wx.NewId(),
                                "K Color Contrast",
                                display_K_menu)
        #subsubmenu Na display
        display_Na_menu = wx.Menu()
        display_Na_menu.AppendRadioItem(self.ID_Na_ContrastHigh, 'high', 'high contrast')
        display_Na_menu.AppendRadioItem(self.ID_Na_ContrastNormal, 'normal', 'low contrast')
        display_Na_menu.AppendRadioItem(self.ID_Na_ContrastLow, 'low', 'high contrast')
        display_Na_menu.AppendRadioItem(self.ID_Na_ContrastUltra, 'ultra', 'ultra high contrast')
        display_menu.AppendMenu(wx.NewId(),
                                "Na Color Contrast",
                                display_Na_menu)

        ##set defaults
        display_K_menu.Check(self.ID_K_ContrastNormal, True) #TODO
        display_Na_menu.Check(self.ID_Na_ContrastNormal, True) #TODO

        self.frame.Bind(wx.EVT_MENU_RANGE, self.OnMenuContrast,
                        id=self.ID_K_ContrastHigh,
                        id2=self.ID_Na_ContrastLow)

        display_menu.AppendSeparator()
        display_menu.Append(self.ID_MarkerA,
                            'Main marker visible',
                            '',
                            wx.ITEM_CHECK)
        display_menu.Append(self.ID_MarkerB,
                            '2nd marker visible',
                            '',
                            wx.ITEM_CHECK)
        display_menu.Append(self.ID_MarkerC,
                            '3rd marker visible',
                            '',
                            wx.ITEM_CHECK)
        self.frame.Bind(wx.EVT_MENU_RANGE,
                        self.OnMenuMarker,
                        id=self.ID_MarkerA,
                        id2=self.ID_MarkerB)
        self.frame.Bind(wx.EVT_MENU_RANGE,
                        self.OnMenuMarker,
                        id=self.ID_MarkerC)

        #submenu Fit
        fit_menu = wx.Menu()
        fit_menu.Append(self.ID_Reload,
                        "Reload image\tF5",
                        "Reload image from file",
                        )
        fit_menu.AppendSeparator()
        sc = fit_menu.AppendCheckItem(self.ID_FitShowContours,
                        'Show Contours',
                        'Show Contour lines')
        sc.Check()
        fit_menu.AppendSeparator()
        
        #subsubmenu Na,1 peak
        fit_Na_menu = wx.Menu()
        fit_Na_menu.AppendRadioItem(
            self.ID_FitNaNone,
            "Na no fit",
            "perform no fit for Sodium",
            )
        fit_Na_menu.AppendRadioItem(
            self.ID_NaImageIntegration,
            "Na image integration",
            "intergrate the raw image",
            )
        fit_Na_menu.AppendRadioItem(
            self.ID_FitNaGauss1D,
            "Na Gauss fit 1D",
            "perform 2 Gauss 1d fits for Sodium",
            )
        fit_Na_menu.AppendRadioItem(
            self.ID_FitNaLorentzGauss1D,
            "Na Lorentz(x)+Gauss(y) fit 1D",
            "perform two 1d fits for Sodium",
            )
        fit_Na_menu.AppendRadioItem(
            self.ID_FitNaGaussGauss1D,
            "Na GaussGauss(x)+Gauss(y) fit 1D",
            "perform 1d fits for Sodium, with 2 gaussian along x",
            )
        fit_Na_menu.AppendRadioItem(
            self.ID_FitNaGauss,
            "Na Gauss fit",
            "perform 2d Gauss fit for Sodium",
            )
        fit_Na_menu.AppendRadioItem(
            self.ID_FitNaGaussSym,
            "Na symmetric Gauss fit",
            "perform symmetric (sx = sy) 2d Gauss fit for Sodium",
            )
        fit_Na_menu.AppendRadioItem(
            self.ID_FitNaGaussBose,
            "Na Bose enhanced Gauss fit",
            "perform 2d Bose enhanced Gauss fit for Sodium",
            )

        fit_Na_menu.AppendRadioItem(
            self.ID_FitNaTF,
            "Na Thomas-Fermi fit",
            "perform 2d Thomas-Fermi fit distribution for Sodium",
            )
        fit_Na_menu.AppendRadioItem(
            self.ID_FitNaTFSlice,
            "Na Thomas-Fermi fit on a Slice",
            "perform 2d Thomas-Fermi fit distribution for Sodium",
            )
        fit_Na_menu.AppendRadioItem(
            self.ID_FitNaBimodal,
            "Na bimodal fit",
            "perform 2d fit of Gaussian + Thomas-Fermi distribution for Sodium",
            )
        fit_Na_menu.AppendRadioItem(
            self.ID_FitNaBimodalSplit,
            "Na bimodal fit - split mxy",
            "perform integrated 1d fit of Gaussian + TF distribution for Sodium, with separate centers",
            )
        fit_Na_menu.AppendRadioItem(
            self.ID_FitNaThomasFermiHole,
            "Na Thomas-Fermi with gaussian hole",
            "perform 2d fit Thomas-Fermi distribution with hole for Sodium",
            )
        fit_Na_menu.AppendRadioItem(
            self.ID_FitNaBimodalIntegrate1dSplit,
            "Na bimodal integrate 1d fit - split mxy",
            "perform 2d fit of Gaussian + Thomas-Fermi distribution for Sodium, with separate centers",
            )
        fit_Na_menu.AppendRadioItem(
            self.ID_FitNaGaussIntegrate1dPinning,
            "Na Gauss fit 1d with hole",
            "perform integrated 1d fit of Gaussian with gaussian hole",
            )
        fit_Na_menu.AppendRadioItem(
            self.ID_FitNaBoseBimodal,
            "Na bose enhanced bimodal fit",
            "perform 2d fit of Bose enhanced Gaussian + Thomas-Fermi distribution for Sodium",
            )
        fit_Na_menu.AppendRadioItem(
            self.ID_FitNaBimodalGaussGauss,
            "Na bimodal Gauss + Gauss fit",
            "perform 2d fit of Gaussian + Gaussian distribution for Sodium",
            )
        fit_Na_menu.AppendRadioItem(
            self.ID_FitNaNGauss,
            "Na Gauss xN fit",
            "perform 2d fit with sum of 3/5 Gaussian distribtions for Sodium",
            )
        fit_menu.AppendMenu(wx.NewId(),
                                "Na",
                                fit_Na_menu)
        
        fit_menu.Append(self.ID_MarkerGuessNa,
                        'Na guess from markers',
                        'Guess peak positions from markers',
                        wx.ITEM_CHECK)
        fit_menu.Append(self.ID_CompensateOffsetNa,
                        'Na compensate background',
                        'Compensate offset of image background',
                        wx.ITEM_CHECK)

        fit_menu.AppendSeparator()
        
        # K
        fit_K_menu = wx.Menu()
        fit_K_menu.AppendRadioItem(
            self.ID_FitKNone,
            "K no fit",
            "perform no fit for Sodium",
            )
        fit_K_menu.AppendRadioItem(
            self.ID_KImageIntegration,
            "K image integration",
            "intergrate the raw image",
            )
        fit_K_menu.AppendRadioItem(
            self.ID_FitKGauss1D,
            "K Gauss fit 1D",
            "perform 2 Gauss 1d fits for Potassium",
            )
        fit_K_menu.AppendRadioItem(
            self.ID_FitKLorentzGauss1D,
            "K Lorentz(x)+Gauss(y) fit 1D",
            "perform two 1d fits for Potassium",
            )
        fit_K_menu.AppendRadioItem(
            self.ID_FitKGaussGauss1D,
            "K GaussGauss(x)+Gauss(y) fit 1D",
            "perform 1d fits for Potassium, with 2 gaussian along x",
            )
        fit_K_menu.AppendRadioItem(
            self.ID_FitKGauss,
            "K Gauss fit",
            "perform 2d Gauss fit for Potassium",
            )
        fit_K_menu.AppendRadioItem(
            self.ID_FitKGaussSym,
            "K symmetric Gauss fit",
            "perform symmetric (sx = sy) 2d Gauss fit for Potassium",
            )
        fit_K_menu.AppendRadioItem(
            self.ID_FitKGaussBose,
            "K Bose enhanced Gauss fit",
            "perform 2d Bose enhanced Gauss fit for Potassium",
            )
        fit_K_menu.AppendRadioItem(
            self.ID_FitKTF,
            "K Thomas-Fermi fit",
            "perform 2d Thomas-Fermi fit distribution for Potassium",
            )
        fit_K_menu.AppendRadioItem(
            self.ID_FitKTFSlice,
            "K Thomas-Fermi fit on a Slice",
            "perform 2d Thomas-Fermi fit distribution for Potassium",
            )
        fit_K_menu.AppendRadioItem(
            self.ID_FitKBimodal,
            "K bimodal fit",
            "perform 2d fit of Gaussian + Thomas-Fermi distribution for Potassium",
            )
        fit_K_menu.AppendRadioItem(
            self.ID_FitKBimodalSplit,
            "K bimodal fit - split mxy",
            "perform 2d fit of Gaussian + Thomas-Fermi distribution for Potassium, with separate centers",
            )
        fit_K_menu.AppendRadioItem(
            self.ID_FitKThomasFermiHole,
            "K Thomas-Fermi with gaussian hole",
            "perform 2d fit Thomas-Fermi distribution with hole for Potassium",
            )
        fit_K_menu.AppendRadioItem(
            self.ID_FitKBimodalIntegrate1dSplit,
            "K bimodal integrate 1d fit - split mxy",
            "perform integrated 1d fit of Gaussian + TF distribution for Potassium, with separate centers",
            )
        fit_K_menu.AppendRadioItem(
            self.ID_FitKGaussIntegrate1dPinning,
            "K Gauss fit 1d with hole",
            "perform integrated 1d fit of Gaussian with gaussian hole",
            )
        fit_K_menu.AppendRadioItem(
            self.ID_FitKBoseBimodal,
            "K bose enhanced bimodal fit",
            "perform 2d fit of Bose enhanced Gaussian + Thomas-Fermi distribution for Sodium",
            )
        fit_K_menu.AppendRadioItem(
            self.ID_FitKBimodalGaussGauss,
            "K bimodal Gauss + Gauss fit",
            "perform 2d fit of Gaussian + Gaussian distribution for Potassium",
            )

        fit_menu.AppendMenu(wx.NewId(),
                                "K",
                                fit_K_menu)
        
        fit_menu.Append(self.ID_MarkerGuessK,
                        'K guess from markers',
                        'Guess peak positions from markers',
                        wx.ITEM_CHECK)
        fit_menu.Append(self.ID_CompensateOffsetK,
                        'K compensate background',
                        'Compensate offset of image background',
                        wx.ITEM_CHECK)
        
        fit_menu.Check(self.ID_FitNaNone, True)
        fit_menu.Check(self.ID_FitKNone, True)
        fit_menu.AppendSeparator()
    
        self.frame.Bind(wx.EVT_MENU,
                        self.OnMenuReload,
                        id=self.ID_Reload)
        self.frame.Bind(wx.EVT_MENU,
                        self.OnMenuFit,
                        id=self.ID_FitShowContours)
        self.frame.Bind(wx.EVT_MENU_RANGE,
                        self.OnMenuFit,
                        id=self.ID_FitNaNone, #first item in fit Na
                        id2=self.ID_FitKBoseBimodal) #last item in fit K
        self.frame.Bind(wx.EVT_MENU_RANGE,
                        self.OnMenuCompensate,
                        id=self.ID_CompensateOffsetNa,
                        id2=self.ID_CompensateOffsetK)
        self.frame.Bind(wx.EVT_MENU_RANGE,
                        self.OnMenuMarkerGuess,
                        id=self.ID_MarkerGuessNa,
                        id2=self.ID_MarkerGuessK)
        
        #submenu: Result
        result_menu = wx.Menu()
        result_menu.Append(self.ID_ResultsPlot,
                           "Plot Results",
                           "Plot Results.")
        
        result_menu.Append(self.ID_ResultsSave,
                           "Save", "Save fit results to file")
        result_menu.Append(self.ID_ResultsSaveAs,
                           "Save as ...", "Save fit results to file")
        result_menu.Append(self.ID_ResultsLoad,
                           "Load...", "Load results from file")

        result_menu.Append(self.ID_ResultsSaveTemplate,
                           "Save Template...", "Save settings to template file")
        result_menu.Append(self.ID_ResultsApplyTemplate,
                           "Apply Template", "Apply settings stored in template file")
        result_menu.Append(self.ID_ResultsNew,
                            "New results from template", "")
        
        
        self.frame.Bind(wx.EVT_MENU, self.OnMenuPlotResults, id=self.ID_ResultsPlot)
        self.frame.Bind(wx.EVT_MENU, self.OnMenuResultsSave, id=self.ID_ResultsSave)
        self.frame.Bind(wx.EVT_MENU, self.OnMenuResultsSave, id=self.ID_ResultsSaveAs)
        self.frame.Bind(wx.EVT_MENU, self.OnMenuResultsLoad, id=self.ID_ResultsLoad)

        self.frame.Bind(wx.EVT_MENU, self.OnMenuResultsSaveTemplate, id=self.ID_ResultsSaveTemplate)
        self.frame.Bind(wx.EVT_MENU, self.OnMenuResultsApplyTemplate, id=self.ID_ResultsApplyTemplate)
        self.frame.Bind(wx.EVT_MENU, self.OnMenuResultsNew, id=self.ID_ResultsNew)
        
        
        #menu: Autosave
        autosave_menu = wx.Menu()
        autosave_menu.Append(self.ID_AutoSaveAbs, 'Save abs image','',wx.ITEM_CHECK)
        autosave_menu.Append(self.ID_AutoSaveRaw, 'Save raw images','',wx.ITEM_CHECK)
        self.frame.Bind(wx.EVT_MENU, self.OnMenuAutoSave, id=self.ID_AutoSaveAbs)
        self.frame.Bind(wx.EVT_MENU, self.OnMenuAutoSave, id=self.ID_AutoSaveRaw)
        
        #menu: Settings
        settings_menu = wx.Menu()
        settings_menu.Append(self.ID_SettingsSave, "Save", "Save Settings")
        settings_menu.Append(self.ID_SettingsLoad, "Load", "Load Settings")

        self.frame.Bind(wx.EVT_MENU, self.OnMenuSettingsSave, id=self.ID_SettingsSave)
        self.frame.Bind(wx.EVT_MENU, self.OnMenuSettingsLoad, id=self.ID_SettingsLoad)
        
        
        #Menus: finalize
        self.menu = wx.MenuBar()

        self.menu.Append(view_menu, "Show")
        self.menu.Append(fit_menu, "Fit")
        self.menu.Append(display_menu, "Display")
        self.menu.Append(result_menu, "Results")
        self.menu.Append(autosave_menu, "AutoSave")
        self.menu.Append(settings_menu, "Settings")

        self.frame.SetMenuBar(self.menu)



        ## create center panel
        self.gridpanel = FitResultTableGrid.TabbedGridPanel(self.frame)
        
        self.Bind(wx.aui.EVT_AUINOTEBOOK_PAGE_CLOSE, self.OnResultsPageClose, self.gridpanel.notebook)
        self.Bind(wx.aui.EVT_AUINOTEBOOK_PAGE_CLOSED, self.OnResultsPageClosed, self.gridpanel.notebook)
        
        self.results_filename = None

        self.mgr.AddPane(self.gridpanel,
                         wx.aui.AuiPaneInfo().
                         Name("results").Caption("Results")
                         .CenterPane().BestSize(wx.Size(800, 100))
                         )

        ## create and initializer image panels

        #imgNa, imgNa = loadimg(self.imagefilename)   #####added 14-12-2012
        #imgK, imgK = loadimg(self.imagefilename.replace(".sis","_1.sis"))   #####added 14-12-2012
        imgK, imgNa = loadimg(self.imagefilename)   #####original 14-12-2012
        imgK = ma.array(imgK, mask= ~ numpy.isfinite(imgK))
        imgNa = ma.array(imgNa, mask= ~ numpy.isfinite(imgNa))
        
        self.K = ImgPanel(self.frame, 'K', imgK, self.roiK[2], fitting.NoFit(self.imaging_pars['K']), pylab.get_cmap(self.imaging_pars['K'].palette))
        self.Na = ImgPanel(self.frame, 'Na', imgNa, self.roiNa[2], fitting.NoFit(self.imaging_pars['Na']), pylab.get_cmap(self.imaging_pars['Na'].palette))

        #TODO: better: reduce constructor, call later show_image
        
        self.mgr.AddPane(self.K, wx.aui.
                         AuiPaneInfo().Name("K").Caption("Potassium").
                         Top().Position(1).BestSize(wx.Size(640, 480))
                         )
        self.mgr.AddPane(self.Na, wx.aui.
                         AuiPaneInfo().Name("Na").Caption("Sodium").
                         Top().Position(2)
                         )
        
        ## create status bar
        
        self.statusbar = StatusBarWx(self.frame)
        self.K.toolbar.set_status_bar(self.statusbar)
        self.Na.toolbar.set_status_bar(self.statusbar)
        self.frame.SetStatusBar(self.statusbar)
        
        ## create toolbars top
        self.tb = wx.ToolBar(self.frame, - 1, wx.DefaultPosition,
                             wx.DefaultSize,
                             wx.TB_FLAT | wx.TB_NODIVIDER
                             )
        
        self.tb.SetToolBitmapSize((24, 24)) #TODO: ???

        self.tb2 = wx.ToolBar(self.frame, - 1, wx.DefaultPosition,
                             wx.DefaultSize,
                             wx.TB_FLAT | wx.TB_NODIVIDER
                             )
        
        self.tb2.SetToolBitmapSize((24, 24)) #TODO: ???
        
        self.tb3 = wx.ToolBar(self.frame, - 1, wx.DefaultPosition,
                             wx.DefaultSize,
                             wx.TB_FLAT | wx.TB_NODIVIDER
                             )
        
        self.tb3.SetToolBitmapSize((24, 24)) #TODO: ???
        
        self.tb4 = wx.ToolBar(self.frame, - 1, wx.DefaultPosition,
                             wx.DefaultSize,
                             wx.TB_FLAT | wx.TB_NODIVIDER
                             )
        
        self.tb4.SetToolBitmapSize((24, 24)) #TODO: ???
        
        self.tb5 = wx.ToolBar(self.frame, - 1, wx.DefaultPosition,
                             wx.DefaultSize,
                             wx.TB_FLAT | wx.TB_NODIVIDER
                             )
        
        self.tb5.SetToolBitmapSize((24, 24)) #TODO: ???
        
        # measurement
        #self.measurement = wx.TextCtrl(self.tb, -1, "", size = (100,-1))
        self.measurement_combobox = wx.ComboBox(self.tb, - 1, "", size=(100, - 1), style=wx.TE_PROCESS_ENTER)
        self.tb.AddControl(self.measurement_combobox)

        self.frame.Bind(wx.EVT_TEXT_ENTER, self.OnChangeMeasurementName, self.measurement_combobox)
        self.frame.Bind(wx.EVT_COMBOBOX, self.OnChangeMeasurementName, self.measurement_combobox)
        self.tb.AddSeparator()

        # Save Image
        self.savebutton = wx.Button(self.tb, self.ID_SaveImage,
                                    "Save Image", size=(100, 30)
                                    )
        self.savebutton.SetBackgroundStyle(wx.BG_STYLE_COLOUR)
        self.tb.AddControl(self.savebutton)
        self.frame.Bind(wx.EVT_BUTTON, self.OnSaveImage, id=self.ID_SaveImage)

        ## Automatic Save Image
        #self.autosave_control = wx.CheckBox(self.tb, self.ID_Autosave,
        #                            "Auto Save Image")
        #self.tb.AddControl(self.autosave_control)
        #self.frame.Bind(wx.EVT_CHECKBOX,
        #                self.OnAutosaveToggle,
        #                id=self.ID_Autosave)
        #self.autosave = False

        # Reload Button
        self.reloadbutton = wx.Button(self.tb, self.ID_ReloadButton,
                                  "Reload", size=(100, 30)
                                  )
        self.tb.AddControl(self.reloadbutton)
        self.frame.Bind(wx.EVT_BUTTON, self.OnReloadButton, id=self.ID_ReloadButton)

        # Automatic reload
        self.autoreload = wx.CheckBox(self.tb, self.ID_Autoreload,
                                      "Auto Reload")
        self.tb.AddControl(self.autoreload)
        self.frame.Bind(wx.EVT_CHECKBOX,
                        self.OnAutoreloadClicked,
                        id=self.ID_Autoreload)

        # Record data
        self.record_data_button = wx.Button(self.tb, self.ID_RecordData,
                                     "Record Data", size=(100, 30)
                                     )
        self.record_data = True
        self.record_data_button.SetBackgroundStyle(wx.BG_STYLE_COLOUR)
        self.record_data_button.SetBackgroundColour(wx.NamedColour('GREEN'))
        self.tb.AddControl(self.record_data_button)
        self.frame.Bind(wx.EVT_BUTTON,
                        self.OnRecordDataButtonClicked,
                        id=self.ID_RecordData)

        ##save image filename display
        #self.tb.AddSeparator()
        #self.saved_image_name = wx.StaticText(self.tb,
        #                                      -1,
        #                                      'image not saved',
        #                                      size = (200, -1),
        #                                      style = wx.ALIGN_LEFT)
        #self.tb.AddControl(self.saved_image_name)

        # Imaging Parameters:
        # Imaging direction: Choicebox
        self.tb2.AddControl(wx.StaticText(self.tb2,
                                          -1,
                                         'Imaging: '))
        choices = []
        for p in self.imaging_parlist:
            choices.append(p['K'].description+" | "+p['Na'].description)
        
        self.imaging_choice = wx.Choice(self.tb2,
                                        self.ID_ImagingChoice,
                                        size=(120, - 1),
                                        choices=choices)
        #self.imaging_choice.Select(0)
        #self.select_imaging(0)
        self.tb2.AddControl(self.imaging_choice)
        self.Bind(wx.EVT_CHOICE, self.OnChoiceImaging, self.imaging_choice)



        choices=["cam0-cam1","cam0-cam0","cam1-cam1"]
        self.imaging_AB =wx.Choice(self.tb2,self.ID_Imaging_AB,size=(80,-1),choices=choices)
        self.tb2.AddControl(self.imaging_AB)
        self.Bind(wx.EVT_CHOICE,self.OnChoice_AB,self.imaging_AB)



        # Expansion times
        self.tb2.AddControl(wx.StaticText(self.tb2,
                                          -1,
                                          " Expansion time K: "))
        self.expansion_time_K = wx.SpinCtrl(self.tb2,
                                            self.ID_ExpansionTimeK,
                                            min=0,
                                            max=500,
                                            initial=10,
                                            name="expansion time K",
                                            size=(50, - 1))
        self.tb2.AddControl(self.expansion_time_K)

        self.tb2.AddControl(wx.StaticText(self.tb2,
                                          -1,
                                          " Na: "))
        self.expansion_time_Na = wx.SpinCtrl(self.tb2,
                                             self.ID_ExpansionTimeNa,
                                             min=0,
                                             max=500,
                                             initial=10,
                                             name="expansion time Na",
                                             size=(50, - 1))
        self.tb2.AddControl(self.expansion_time_Na)

        self.Bind(wx.EVT_SPINCTRL, self.OnChangeExpansionTime, self.expansion_time_K)
        self.Bind(wx.EVT_TEXT_ENTER, self.OnChangeExpansionTime, self.expansion_time_K)
        self.Bind(wx.EVT_SPINCTRL, self.OnChangeExpansionTime, self.expansion_time_Na)
        self.Bind(wx.EVT_TEXT_ENTER, self.OnChangeExpansionTime, self.expansion_time_Na)                                           
        #maximum optical density entries
        self.tb2.AddControl(wx.StaticText(self.tb2,
                                          -1,
                                          " Max. OD K: "))
        self.entry_optical_density_max_K = wx.TextCtrl(self.tb2,
                                                 self.ID_OpticalDensityMaxK,
                                                 '0.0',
                                                 size=(50, - 1),
                                                 name="max. OD K",
                                                 style=wx.TE_PROCESS_ENTER,
                                                 )
        self.tb2.AddControl(self.entry_optical_density_max_K)
        
        self.tb2.AddControl(wx.StaticText(self.tb2,
                                          -1,
                                          " Na: "))
        self.entry_optical_density_max_Na = wx.TextCtrl(self.tb2,
                                                 self.ID_OpticalDensityMaxNa,
                                                 '0.0',
                                                 size=(50, - 1),
                                                 name="max. OD Na",
                                                 style=wx.TE_PROCESS_ENTER,
                                                 )
        self.tb2.AddControl(self.entry_optical_density_max_Na)
        
        self.Bind(wx.EVT_TEXT_ENTER, self.OnChangeMaximumOpticalDensity, self.entry_optical_density_max_Na)
        self.Bind(wx.EVT_TEXT_ENTER, self.OnChangeMaximumOpticalDensity, self.entry_optical_density_max_K)
                


        self.tb2.AddControl(wx.StaticText(self.tb2,
                                          -1,
                                          " Rot K: "))
        self.rotation_K = wx.TextCtrl(self.tb2,
                                                 self.ID_RotateK,
                                                 '0.0',
                                                 size=(50, - 1),
                                                 name="Rotation K",
                                                 style=wx.TE_PROCESS_ENTER,
                                                 )
        self.tb2.AddControl(self.rotation_K)
        
        self.tb2.AddControl(wx.StaticText(self.tb2,
                                          -1,
                                          " Na: "))
        self.rotation_Na = wx.TextCtrl(self.tb2,
                                                 self.ID_RotateNa,
                                                 '0.0',
                                                 size=(50, - 1),
                                                 name="Rotation Na",
                                                 style=wx.TE_PROCESS_ENTER,
                                                 )
        self.tb2.AddControl(self.rotation_Na)
        
        self.Bind(wx.EVT_TEXT_ENTER, self.OnChangeRotation, self.rotation_Na)
        self.Bind(wx.EVT_TEXT_ENTER, self.OnChangeRotation, self.rotation_K)





        self.tb5.AddControl(wx.StaticText(self.tb5,
                                          -1,
                                          " Vmin K: "))
        self.vmin_K = wx.TextCtrl(self.tb5,
                                                 self.ID_VminK,
                                                 '-0.05',
                                                 size=(50, - 1),
                                                 name="Vmin K",
                                                 style=wx.TE_PROCESS_ENTER,
                                                 )
        self.tb5.AddControl(self.vmin_K)
        
        self.tb5.AddControl(wx.StaticText(self.tb5,
                                          -1,
                                          " Na: "))
        self.vmin_Na = wx.TextCtrl(self.tb5,
                                                 self.ID_VminNa,
                                                 '-0.05',
                                                 size=(50, - 1),
                                                 name="Vmin Na",
                                                 style=wx.TE_PROCESS_ENTER,
                                                 )
        self.tb5.AddControl(self.vmin_Na)
        
        self.Bind(wx.EVT_TEXT_ENTER, self.OnChangeVmin, self.vmin_Na)
        self.Bind(wx.EVT_TEXT_ENTER, self.OnChangeVmin, self.vmin_K)
        
        
        self.tb5.AddControl(wx.StaticText(self.tb5,
                                          -1,
                                          " Vmax K: "))
        self.vmax_K = wx.TextCtrl(self.tb5,
                                                 self.ID_VmaxK,
                                                 '1.5',
                                                 size=(50, - 1),
                                                 name="Vmax K",
                                                 style=wx.TE_PROCESS_ENTER,
                                                 )
        self.tb5.AddControl(self.vmax_K)
        
        self.tb5.AddControl(wx.StaticText(self.tb5,
                                          -1,
                                          " Na: "))
        self.vmax_Na = wx.TextCtrl(self.tb5,
                                                 self.ID_VmaxNa,
                                                 '1.5',
                                                 size=(50, - 1),
                                                 name="Vmax Na",
                                                 style=wx.TE_PROCESS_ENTER,
                                                 )
        self.tb5.AddControl(self.vmax_Na)
        
        self.Bind(wx.EVT_TEXT_ENTER, self.OnChangeVmax, self.vmax_Na)
        self.Bind(wx.EVT_TEXT_ENTER, self.OnChangeVmax, self.vmax_K)



        self.tb2.AddSeparator()
        self.autoscaleclimbutton = wx.Button(self.tb2, self.ID_AutoscaleClim,
                                    "Auto Vmin/max", size=(100, 30)
                                    )
        self.autoscaleclimbutton.SetBackgroundStyle(wx.BG_STYLE_COLOUR)
        self.tb2.AddControl(self.autoscaleclimbutton)
        self.frame.Bind(wx.EVT_BUTTON, self.OnAutoscaleClim, id=self.ID_AutoscaleClim)






        self.tb3.AddControl(wx.StaticText(self.tb3,
                                          -1,
                                          "K residuals:"))
        self.differenceK = wx.CheckBox(self.tb3, self.ID_DifferenceK,
                                      "")
        
        self.tb3.AddControl(self.differenceK)
        
        
        self.tb3.AddControl(wx.StaticText(self.tb3,
                                          -1,
                                          " filter:"))
        self.filterK = wx.CheckBox(self.tb3, self.ID_FilterK,
                                      "")
        self.tb3.AddControl(self.filterK)
        
        
        
        
        
        self.tb3.AddControl(wx.StaticText(self.tb3,
                                          -1,
                                          " filter sigma: "))
        self.entry_filter_sigma_K = wx.TextCtrl(self.tb3,
                                                 self.ID_FilterSigmaK,
                                                 str(self.K.filter_sigma),
                                                 size=(50, - 1),
                                                 name="sigma K",
                                                 style=wx.TE_PROCESS_ENTER,
                                                 )
        self.tb3.AddControl(self.entry_filter_sigma_K)
        
        
        self.tb3.AddControl(wx.StaticText(self.tb3,
                                          -1,
                                          " scale: "))
        self.entry_filter_scalar_K = wx.TextCtrl(self.tb3,
                                                 self.ID_FilterScalarK,
                                                 str(self.K.filter_scalar),
                                                 size=(50, - 1),
                                                 name="scale K",
                                                 style=wx.TE_PROCESS_ENTER,
                                                 )
        self.tb3.AddControl(self.entry_filter_scalar_K)
        
        self.tb3.AddControl(wx.StaticText(self.tb3,
                                          -1,
                                          " reference: "))
        self.referenceK = wx.CheckBox(self.tb3, self.ID_ReferenceK,
                                      "")
        self.tb3.AddControl(self.referenceK)
        
        
        
        self.tb3.AddSeparator()
        self.tb3.AddControl(wx.StaticText(self.tb3,
                                          -1,
                                          " Na residuals: "))
        self.differenceNa = wx.CheckBox(self.tb3, self.ID_DifferenceNa,
                                      "")

        self.tb3.AddControl(self.differenceNa)
        
        self.tb3.AddControl(wx.StaticText(self.tb3,
                                          -1,
                                          " filter: "))
        self.filterNa = wx.CheckBox(self.tb3, self.ID_FilterNa,
                                      "")
        self.tb3.AddControl(self.filterNa)
        
        self.tb3.AddControl(wx.StaticText(self.tb3,
                                          -1,
                                          " filter sigma: "))
        self.entry_filter_sigma_Na = wx.TextCtrl(self.tb3,
                                                 self.ID_FilterSigmaNa,
                                                 str(self.Na.filter_sigma),
                                                 size=(50, - 1),
                                                 name="sigma Na",
                                                 style=wx.TE_PROCESS_ENTER,
                                                 )
        self.tb3.AddControl(self.entry_filter_sigma_Na)
        
        
        self.tb3.AddControl(wx.StaticText(self.tb3,
                                          -1,
                                          " scale: "))
        self.entry_filter_scalar_Na = wx.TextCtrl(self.tb3,
                                                 self.ID_FilterScalarNa,
                                                 str(self.Na.filter_scalar),
                                                 size=(50, - 1),
                                                 name="scale Na",
                                                 style=wx.TE_PROCESS_ENTER,
                                                 )
        self.tb3.AddControl(self.entry_filter_scalar_Na)
        
        self.tb3.AddControl(wx.StaticText(self.tb3,
                                          -1,
                                          " reference: "))
        self.referenceNa = wx.CheckBox(self.tb3, self.ID_ReferenceNa,
                                      "")
        self.tb3.AddControl(self.referenceNa)
        
        
        self.frame.Bind(wx.EVT_CHECKBOX,
                        self.OnDifferenceSelectK,
                        id=self.ID_DifferenceK)
        
        
        self.frame.Bind(wx.EVT_CHECKBOX,
                        self.OnDifferenceSelectNa,
                        id=self.ID_DifferenceNa)

        
        self.frame.Bind(wx.EVT_CHECKBOX,
                        self.OnFilterSelectK,
                        id=self.ID_FilterK)

        self.frame.Bind(wx.EVT_CHECKBOX,
                        self.OnFilterSelectNa,
                        id=self.ID_FilterNa)


        self.frame.Bind(wx.EVT_CHECKBOX,
                        self.OnReferenceSelect,
                        id=self.ID_ReferenceK)


        self.frame.Bind(wx.EVT_CHECKBOX,
                        self.OnReferenceSelect,
                        id=self.ID_ReferenceNa)

        

        
        self.Bind(wx.EVT_TEXT_ENTER, self.OnChangeFilterSigma, self.entry_filter_sigma_Na)
        self.Bind(wx.EVT_TEXT_ENTER, self.OnChangeFilterSigma, self.entry_filter_sigma_K)
        
        self.Bind(wx.EVT_TEXT_ENTER, self.OnChangeFilterScalar, self.entry_filter_scalar_Na)
        self.Bind(wx.EVT_TEXT_ENTER, self.OnChangeFilterScalar, self.entry_filter_scalar_K)
                



        
        self.savebitmapbutton = wx.Button(self.tb4, self.ID_SaveBitmap,
                                    "Save Bitmaps", size=(100, 30)
                                    )
        self.savebitmapbutton.SetBackgroundStyle(wx.BG_STYLE_COLOUR)
        self.tb4.AddControl(self.savebitmapbutton)
        self.frame.Bind(wx.EVT_BUTTON, self.OnSaveBitmap, id=self.ID_SaveBitmap)

        
        self.tb4.AddControl(wx.StaticText(self.tb4,
                                          -1,
                                         ' Colors K: '))
        
        self.maps=[map for map in pylab.cm.datad if not map.endswith("_r")]
        self.maps.sort()
        
        self.cmap_choiceK = wx.Choice(self.tb4,
                                        self.ID_CMapChoiceK,
                                        size=(75, - 1),
                                        choices=self.maps)
        self.tb4.AddControl(self.cmap_choiceK)
        self.Bind(wx.EVT_CHOICE, self.OnChoiceCMap, self.cmap_choiceK)
        
        self.tb4.AddControl(wx.StaticText(self.tb4,
                                          -1,
                                         ' Na: '))
        self.cmap_choiceNa = wx.Choice(self.tb4,
                                        self.ID_CMapChoiceNa,
                                        size=(75, - 1),
                                        choices=self.maps)
        self.tb4.AddControl(self.cmap_choiceNa)
        self.Bind(wx.EVT_CHOICE, self.OnChoiceCMap, self.cmap_choiceNa)

        self.createreference = wx.Button(self.tb4, self.ID_CreateReference,
                                  "Create reference", size=(100, 30)
                                  )
        self.tb4.AddControl(self.createreference)
        self.frame.Bind(wx.EVT_BUTTON, self.OnCreateReference, id=self.ID_CreateReference)



        #finalize toolbars
        self.tb.Realize()
        self.tb2.Realize()
        self.tb3.Realize()
        self.tb4.Realize()
        self.tb5.Realize()
        
        self.mgr.AddPane(self.tb,
                         wx.aui.AuiPaneInfo().
                         Name("toolbar").
                         ToolbarPane().Top().Row(1).Position(1))


        self.mgr.AddPane(self.tb2,
                         wx.aui.AuiPaneInfo().
                         Name("toolbar2").
                         ToolbarPane().Top().Row(1).Position(2))
        self.mgr.AddPane(self.tb3,
                         wx.aui.AuiPaneInfo().
                         Name("toolbar3").
                         ToolbarPane().Top().Row(2).Position(1))
        self.mgr.AddPane(self.tb4,
                         wx.aui.AuiPaneInfo().
                         Name("toolbar4").
                         ToolbarPane().Top().Row(2).Position(2))
                         
        self.mgr.AddPane(self.tb5,
                         wx.aui.AuiPaneInfo().
                         Name("toolbar5").
                         ToolbarPane().Top().Row(2).Position(3))

        #ImageTree
        self.savedimagestree = ImageTree.TreeModelSync(os.path.normpath(settings.imagesavepath)) #TODO
        self.imagetreepanel = ImageTree.ImageTreePanel(self.frame,
                                             self.savedimagestree)
        self.selectedsavedimage = None
        self.mgr.AddPane(self.imagetreepanel,
                         wx.aui.AuiPaneInfo().
                         Name("imagetree").
                         Caption("Saved Images").
                         Right().BestSize(wx.Size(200, 400))
                         )
        self.imagetreepanel.Bind(wx.EVT_TREE_ITEM_ACTIVATED, self.OnActivateImageTree)
        self.imagetreepanel.Bind(wx.EVT_TREE_ITEM_RIGHT_CLICK, self.OnRightClickImageTree)
        

        #ImagePanels for saved images
        self.savedimageK = ImagePanel.CamImageMarkersPanel(self.frame)
        self.savedimageNa = ImagePanel.CamImageMarkersPanel(self.frame)

        for panel in [self.savedimageK, self.savedimageNa]:
            panel.toolbar.AddSeparator()
            panel.toolbar.AddLabelTool(self.ID_AnalyzeImageButton,
                                       'Analyze',
                                       wx.Bitmap(os.path.join(settings.bitmappath,
                                                              'cam24.png'),
                                                 wx.BITMAP_TYPE_PNG),
                                       shortHelp = 'Analyze',
                                       longHelp = 'load image in cam for analyzing')
            panel.toolbar.Realize()
            panel.Bind(wx.EVT_TOOL, self.OnAnalyzeImage, id=self.ID_AnalyzeImageButton)

        self.mgr.AddPane(self.savedimageK,
                         wx.aui.AuiPaneInfo().
                         Name("savedimageK").
                         Caption("saved image of K").
                         Float().Dockable(False).
                         BestSize(wx.Size(400, 300))
                         )
        self.mgr.AddPane(self.savedimageNa,
                         wx.aui.AuiPaneInfo().
                         Name("savedimageNa").
                         Caption("saved image of Na").
                         Float().Dockable(False).
                         BestSize(wx.Size(400, 300))
                         )
        
        ## create some perspectives
        self.mgr.GetPane("K").Show()
        self.mgr.GetPane("Na").Show()
        self.mgr.GetPane("results").Show()
        self.mgr.GetPane("Toolbar").Show()
        self.mgr.GetPane("savedimageNa").Hide()
        self.mgr.GetPane("savedimageK").Hide()
        self.perspective_full = self.mgr.SavePerspective()

        self.mgr.GetPane("K").Hide()
        self.mgr.GetPane("Na").Show().Left()
        self.perspective_Na = self.mgr.SavePerspective()

        self.mgr.GetPane("Na").Hide()
        self.mgr.GetPane("K").Show().Left()
        self.perspective_K = self.mgr.SavePerspective()

        self.mgr.LoadPerspective(self.perspective_full)

        ## Reload Event
        self.frame.Bind(EVT_RELOAD_IMAGE, self.OnReloadEvent)

        ## Save Image Event
        self.frame.Bind(EVT_SAVE_IMAGE, self.OnSaveImageEvent)

        ##Fit update event
        self.frame.Bind(EVT_FIT_RESULTS, self.OnFitResultsEvent)
        
        #tweak scaling
        self.Na.axhprof.set_xlim(self.Na.roi.xmin - 100, self.Na.roi.xmax + 100)
        self.Na.toolbar.push_current()

        ##initialize








        ######default settings 2013-12-03######################################
        fit_menu.Check(self.ID_FitNaGauss, True)
        fit_menu.Check(self.ID_FitKGauss, True)
        fit_menu.Check(self.ID_CompensateOffsetNa, True)
        fit_menu.Check(self.ID_CompensateOffsetK, True)
        
        self.Na.fit = fitting.Gauss2d(self.imaging_pars['Na'])
        self.Na.update()
        self.K.fit = fitting.Gauss2d(self.imaging_pars['K'])
        self.K.update()
        self.Na.set_compensate_offset(True)
        self.K.set_compensate_offset(True)
        self.Na.OnAutoscale(None)
        self.K.OnAutoscale(None)
        
        self.autoreload.SetValue(True)
        self.OnAutoreloadClicked(self.autoreload)

        
        self.plotpanels = dict()
        self.ding = ding.Ding()
        self.InitSettings()

        self.mgr.Update()
        self.frame.Show(True)

        return True


############# added 17-12-2012
    def ConcatenateSis(self):
        #self.imagefilename
        
        img0 = read(self.imagefilename.replace(".sis","_0.sis"))
        img1 = read(self.imagefilename.replace(".sis","_1.sis"))
        h, w = img0.shape

        img=numpy.concatenate((img0[h/2:], img1[h/2:]))

        img.shape=img0.shape

        #img=numpy.concatenate([img0, img1])
                
        write_raw_image(self.imagefilename, img)
##############


    def InitSettings(self):
        """Initialize all values and settings which is not yet possible to do 
        in OnInit since some of the controls do not yet exist"""
        self.Na.show_contours = True
        self.K.show_contours = True
        
        self.cmap_choiceK.SetStringSelection(self.imaging_pars['K'].palette)
        self.cmap_choiceNa.SetStringSelection(self.imaging_pars['Na'].palette)
        self.imaging_choice.Select(0)
        self.select_imaging(0)
        self.imaging_AB.Select(0)

        
        self.measurements = []
        self.select_measurement("results")
        self.sync_imaging_pars_expansion_time()
        imgK, imgNa = loadimg(self.imagefilename)
        ############# added 17-12-2012
        self.ConcatenateSis()
        self.fileconcatenatethread0 = filewatch.FileChangeNotifier(self.imagefilename.replace(".sis","_0.sis"),
                                                    callback=self.ConcatenateSis)
        self.fileconcatenatethread1 = filewatch.FileChangeNotifier(self.imagefilename.replace(".sis","_1.sis"),
                                                    callback=self.ConcatenateSis)
        self.fileconcatenatethread0.start()
        #self.fileconcatenatethread1.start()
        #############
        
        
        
    def UpdateResults(self):
        self.results.UpdateResults({'K': self.K.fitpars, 'Na': self.Na.fitpars})

    def UpdateResultsFilename(self, filename):
        """Update entry filename in results (after saving or reloading from file)"""
        head, tail = os.path.split(filename)
        self.results.UpdateFilename(tail)
        
    def UpdateResultsTimestamp(self, timestamp):
        """Update entry filename in results (after saving or reloading from file)"""
        self.results.UpdateTimestamp(timestamp)

    def OnReloadEvent(self, event):
        """Reload image. load image (default location if
        event.filename == None), store results into event.target"""

        if event.filename is not None:
            full_path = self.savedimagestree.find_file(event.filename)
            print full_path
            if full_path:
                filename = full_path
            else:
                print "can't find image with filename %s!"%event.filename
                return
        else:
            filename = self.imagefilename
            
        self.results.activate() #activate results table on reload
        self.load_image(filename, event.target)
        self.savebutton.SetBackgroundColour(wx.NamedColor("RED"))
        self.savebutton.Refresh() #necessary?
        
        #self.saved_image_name.SetLabel('image not saved') #TODO: ??
        
        self.statusbar.StatusText = 'image not saved'

        if (self.autosave_abs or self.autosave_raw) and event.target is None:
            wx.PostEvent(self.frame, SaveImageEvent())
        event.Skip()

    def load_image(self, filename, target = None):

        """load image, analyze image. Store results in target. Note
        results are created asynchronously in parallel threads,
        received by FitResultEvent."""
        #imgNa, imgNa = loadimg(filename)   #####added 14-12-2012
        #imgK, imgK = loadimg(filename.replace(".sis","_1.sis"))   #####added 14-12-2012
        imgK, imgNa = loadimg(filename)   #####original 14-12-2012

        imgK_old=imgK
        imgNa_old=imgNa
        
        if self.imaging_AB.GetSelection()==0:
            imgK=imgK_old
            imgNa=imgNa_old
        if self.imaging_AB.GetSelection()==1:
            imgK=imgK_old
            imgNa=imgK_old
        if self.imaging_AB.GetSelection()==2:
            imgK=imgNa_old
            imgNa=imgNa_old
        
        #loading new image also means appending row
        if self.record_data and target is None:
            self.results.AppendRows(1) #TODO: this should be done somewhere else

        self.K.analyze_image(imgK, target)
        self.Na.analyze_image(imgNa, target)

        
    def OnFitResultsEvent(self, event):
        species = event.source
        target = event.target
        fitpars = event.fitpars
        if target is None:
            self.results.UpdateResults({species: fitpars})
        else:
            if target['name'] != self.results.name:
                print "target measurement '%s' does not match active measurement '%s'"%(target['name'], self.results.name)
            else:
                self.results.UpdateResults(data = {species: fitpars},
                                           row = target['row'])

        event.Skip()
        self.UpdateResiduals()
        self.results_autosave()



    def OnAutoreloadClicked(self, event):
        if event.IsChecked():
            #activate automatic reload
            self.filewatchthread = filewatch.FileChangeNotifier(self.imagefilename,    ############commented on 01-02-2013
                                                callback=self.CreateReloadEvent)
            self.filewatchthread.start()
        else:
            self.filewatchthread.keeprunning = False
            self.filewatchthread.join()
            

    def OnRecordDataButtonClicked(self, event):
        if self.record_data:
            self.record_data_button.SetBackgroundColour(wx.NamedColor("RED"))
            #self.record_data_button.SetForegroundColour(wx.NamedColor("RED"))
            self.record_data = False
            self.results.record = False
        else:
            self.record_data_button.SetBackgroundColour(wx.NamedColor("GREEN"))
            #self.record_data_button.SetForegroundColour(wx.NamedColor("GREEN"))
            self.record_data = True
            self.results.record = True

        self.record_data_button.Refresh()
        self.gridpanel.Refresh()

    def PostReloadImageEvent(self):
        """post an ReloadImageEvent with empty arguments: load default
        image, append results to active measurement""" 
        wx.PostEvent(self.frame, ReloadImageEvent(filename = None, target = None))

    def AcousticSignal(self):
        pass
        #self.ding.play()

    def CreateReloadEvent(self):
        self.AcousticSignal()
        self.PostReloadImageEvent()

    def OnMenuSettingsSave(self, event):
        import shelve

        S = shelve.open(os.path.join(settings.basedir,
                                            'settings/settings_cam'))
        try:
            #Window layout
            S['perspective'] = self.mgr.SavePerspective()

            #regions of interest
            S['roiK'] = self.roiK
            S['roiNa'] = self.roiNa

            #imaging pars
            S['imaging_selection'] = self.imaging_choice.GetSelection()

            #markers
            markersKpositions = []
            for k, marker in enumerate(self.K.markers):
                markersKpositions.append(marker.position)

            markersNapositions = []
            for k, marker in enumerate(self.Na.markers):
                markersNapositions.append(marker.position)

            S['markersKpositions'] = markersKpositions
            S['markersNapositions'] = markersNapositions

            S['imaging_parlist'] = self.imaging_parlist
            S['imaging_choice'] = self.imaging_choice.Selection
            
                
        finally:
            S.close()
        
    def OnMenuSettingsLoad(self, event):
        import shelve
        #with shelve.open('settings') as settings:
        S = shelve.open(os.path.join(settings.basedir,
                                            'settings/settings_cam'))
        try:
            #perspectives
            self.mgr.LoadPerspective(S['perspective'])

            #region of interests
            self.roiK = S['roiK']
            self.roiNa = S['roiNa']

            #imaging settings
            imaging_sel = S['imaging_selection']
            self.imaging_choice.SetSelection(imaging_sel)
            self.select_imaging(imaging_sel)
            self.select_roi(imaging_sel)

            #markers
            mKpos = S['markersKpositions']
            for k, pos in enumerate(mKpos):
                self.K.markers[k].position = pos

            mNapos = S['markersNapositions']
            for k, pos in enumerate(mNapos):
                self.Na.markers[k].position = pos
                
            self.imaging_parlist = S['imaging_parlist']
            imagingsel = S['imaging_choice'] 
            
            self.imaging_choice.Selection = imagingsel
            self.select_roi(imagingsel)
            self.imaging_parlist = S['imaging_parlist']
            
            self.select_imaging(imagingsel)
             
        finally:
            S.close()
        
        self.K.center_roi()
        self.Na.center_roi()
            
        self.K.update()
        self.Na.update()

    

    #def reload(self, mplevent):
    #    imgK, imgNa = loadimg(self.filename)
    #    if mplevent.inaxes == self.Na.axtop:
    #        self.Na.show_image(imgNa)
    #        self.K.show_image(imgK)
    #    elif mplevent.inaxes == self.K.axtop:
    #        self.K.show_image(imgK)
    #        self.Na.show_image(imgNa)
    #    else:
    #        self.K.show_image(imgK)
    #        self.Na.show_image(imgNa)
    #    #TODO: share following code with OnReloadEvent ?
    #    self.results.AppendRows()
    #    self.UpdateResults()
    #    self.savebutton.SetBackgroundColour(wx.NamedColor("RED"))


    def OnChangeMeasurementName(self, event):
        name = event.GetString()
        self.select_measurement(name)

    def new_measurement(self, name):
        """create notebook panel, add name to
         combo box list and set combo box
         value
        """
        self.gridpanel.addpage(name)
        self.measurement_combobox.Append(name)
        self.measurement_combobox.Value = name
        self.measurements.append(name)
        
    def delete_measurement(self, name):
        #remove from combobox
        index = self.measurements.index(name)
        self.measurement_combobox.Delete(index)
        #remove from internal list
        self.measurements.remove(name)

        #remove plotpanel
        panelname = name + ' plot' #see plot_results
        plotpane = self.mgr.GetPane(panelname)
        if plotpane.IsOk():
            datapanel = self.plotpanels.pop(panelname)
            self.mgr.DetachPane(datapanel)
            #del datapanel
        
    def select_measurement(self, name):
        """select measurement. if measurement does not exist, add new
        page, update control.
        """
        if self.measurement_combobox.FindString(name) == wx.NOT_FOUND:
            self.new_measurement(name)
            
        index = self.measurements.index(name)
        self.measurement_combobox.SetSelection(index)
        self.gridpanel.notebook.SetSelection(index)

        #deactivate Table
        try:
            self.results.activate(False) #will fail if self.results is
                                         #not yet initialized
        except AttributeError:
            pass 

        #change short reference 'results' and 
        self.results = self.gridpanel.pages[index].Table
        self.results.name = name
        #self.results.activate() #activation done on first reload

        #mark active page, first remove mark
        for i in range(self.gridpanel.notebook.GetPageCount()):
            self.gridpanel.notebook.SetPageBitmap(i, wx.NullBitmap)
        self.gridpanel.notebook.SetPageBitmap(index,
                                              wx.Bitmap(
                                                  os.path.join(settings.bitmappath,
                                                               'star-small.png'),
                                                  wx.BITMAP_TYPE_PNG)
                                              )
        self.gridpanel.Refresh()

    def create_plotpanel(self, measurement_name, panelname):
        """create plotpanel, add to AUI Manager, establish connections
        for automatic updates, store it."""

        plotpanel = DataPanel.DataPlotPanel(self.frame, measurement_name)
        self.mgr.AddPane(plotpanel,
                         wx.aui.AuiPaneInfo()
                         .Name(panelname)
                         .Caption(measurement_name)
                         .BestSize(wx.Size(400, 300))
                         .Float()
                         .Dockable(False)
                         )
        self.mgr.Update()

        #plotpanel.set_table(self.results)
        self.results.add_observer(plotpanel)
        return plotpanel

    def OnMenuPlotResults(self, event):
        """Handle selection of 'plot results' menu. Calls L{create_plotpanel}
        """
        self.plot_results()

    def plot_results(self):
        """Create plot panel
        if necessary, show if hidden, update plotpanel.
        """

        measurement_name = self.results.name
        panelname = measurement_name + ' plot'
        
        #try to get plot panel (AUIPaneInfo) associated with current measurement
        plotpane = self.mgr.GetPane(panelname)
        if not plotpane.IsOk(): #plot panel does not yet exist
            self.plotpanels[panelname] = self.create_plotpanel(measurement_name, panelname)
            plotpane = self.mgr.GetPane(panelname)

        if not plotpane.IsShown():    
            plotpane.Show()
            self.mgr.Update()
            self.plotpanels[panelname].refresh()
            
        #TODO: plotpanel skips updates if not visible
        #self.plotpanels[name].update()
        self.plotpanels[panelname].draw()
        

    def OnActivateImageTree(self, event):
        index = self.imagetreepanel.treectrl.GetIndexOfItem(event.Item)
        filename = self.savedimagestree.GetItemFile(index)

        if len(index)<3:
            #not on image
            return

        if self.savedimagestree.check_day(index[0]):
            print "refresh day"
            self.imagetreepanel.treectrl.RefreshItems()

        self.displayedsavedimage = filename

        #imgK, imgK = ImagePanel.loadimg(filename)       ######added 14-12-2012
        #imgNa, imgNa = ImagePanel.loadimg(filename.replace(".sis","_1.sis"))       ######added 14-12-2012
        imgK, imgNa = ImagePanel.loadimg(filename)       ######original 14-12-2012
        self.savedimageK.show_image(imgK, description=filename)
        self.savedimageNa.show_image(imgNa, description=filename)

        paneimgna = self.mgr.GetPane('savedimageNa')
        paneimgk = self.mgr.GetPane('savedimageK')
        if not paneimgna.IsShown() and not paneimgk.IsShown():
            paneimgna.Float().Show()
            paneimgk.Float().Show()
            self.mgr.Update()

        next = self.imagetreepanel.treectrl.GetNextSibling(event.Item)
        if next.IsOk():
            self.imagetreepanel.treectrl.SelectItem(next)

    def OnRightClickImageTree(self, event):
        index = self.imagetreepanel.treectrl.GetIndexOfItem(event.Item)
        if len(index)>=3:
            filename = self.savedimagestree.GetItemFile(index)
            self.selectedsavedimage = filename
        else:
            self.selectedsavedimage = None

        menu = wx.Menu()
        menu.Append(self.ID_ImageTreeRescan, 'Rescan')
        if len(index)>=1:
            menu.Append(-1, "Day Menu")
        if len(index)>=2:
            menu.Append(-1, "Measurement Menu")
        if len(index)>=3:
            menu.Append(self.ID_ImageTreeReloadSavedImage, "Reload Image")
        self.imagetreepanel.Bind(wx.EVT_MENU, self.OnReloadSavedImage, id = self.ID_ImageTreeReloadSavedImage)
        self.imagetreepanel.Bind(wx.EVT_MENU, self.OnRescanImageTree, id = self.ID_ImageTreeRescan)
        self.imagetreepanel.PopupMenu(menu)
        menu.Destroy()

    def OnAnalyzeImage(self, event):
        """Load and analyze images. Bound to button in saved image panels"""
        self.load_image(self.displayedsavedimage)
        self.UpdateResultsFilename(self.displayedsavedimage)
        self.statusbar.StatusText = 'reload: ' + self.displayedsavedimage

    def OnReloadSavedImage(self, event):
        """Load and analyze images. Bound to popup menu in image tree"""
        self.load_image(self.selectedsavedimage)
        self.UpdateResultsFilename(self.selectedsavedimage)
        self.statusbar.StatusText = 'reload: ' + self.selectedsavedimage

    def OnRescanImageTree(self, event):
        print "rescan files", self.imagetreepanel.treemodel.root
        self.imagetreepanel.treemodel.createfiletree()
        self.savedimagestree = self.imagetreepanel.treemodel
        self.imagetreepanel.treectrl.RefreshItems()
        
    def OnResultsPageClose(self, event):
        page = self.gridpanel.notebook.GetPage(event.Selection)
        results = page.Table
        need_save = results.modified
        name = results.name
        
        if need_save:
            dlg = wx.MessageDialog(self.frame,
                               "Data changes are not yet saved.\n"+\
                                "Do you want to save data?\n"+\
                                "(otherwise they are lost, forever!)",
                               caption = "Close data window",
                               style = wx.YES_NO | wx.CANCEL | wx.ICON_EXCLAMATION,
                               )
            answer = dlg.ShowModal()
            dlg.Destroy()
            if answer == wx.ID_YES:
                success = self.results_save(filename = results.filename, results = results, forcesave = True)
                if success:
                    self.delete_measurement(name)
                else:
                    event.Veto()
            elif answer == wx.ID_NO:
                print "don't save data, it's your choice"
                self.delete_measurement(name)
            else: #CANCEL
                event.Veto()
        else:
            self.delete_measurement(name)
        
    def OnResultsPageClosed(self, event):
        self.select_measurement("results")
        print "OnResultsPageClosed"
        
"""
def test_img_from_disk():

    d = shelve.open('sis.pkl')

    print d.keys()
    img1 = d['img1']
    img2 = d['img2']
    bkg = d['img3']
    d.close()

    img1 = (img1 - bkg).astype(Float)
    img2 = (img2 - bkg).astype(Float)

    den = - log(img1 / img2)

    imgK = den[0:1040, :]
    imgNa = den[1040:2080, :]

    img = ma.array(imgNa, mask= ~ isfinite(imgNa))
    img.set_fill_value(0)

    #img = ma.masked_inside(imgNa, -0.1, 3)


    roi1 = ROI(450, 650, 350, 480)
    roi2 = ROI(450, 650, 350, 481)

    imgdisp = ImgDisp(img, roi1, fitting.Fitting(), fig=1)
    #imgdisp2 = ImgDisp(ma.array(imgK, mask = ~isfinite(imgK)), roi2, Fitting(), fig=2)


    show()
"""

def run_cam():
    gui = ImgAppAui(redirect=False)
    gui.MainLoop()    
    return gui

"""
def profile():
    import cProfile
    cProfile.run('run()', stats)
"""    
if __name__ == '__main__':
    gui = run_cam()
    
