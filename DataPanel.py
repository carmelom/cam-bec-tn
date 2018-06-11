#!/usr/bin/python
#-*- coding: latin-1 -*-
"""Plot measurement results."""

from __future__ import with_statement

import wx
#import math
import pylab, numpy
import matplotlib
matplotlib.pyplot.rcParams['axes.grid'] = True
matplotlib.pyplot.rcParams['image.interpolation'] = 'nearest'
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wx import StatusBarWx
from matplotlib.figure import Figure

import matplotlib.axis

import FitResultTableGrid
reload(FitResultTableGrid)

import toolbar
reload(toolbar)
from toolbar import NavigationToolbar2Wx

from avgdata import avgdata
import os.path
import settings
import time
import math
import types
import toolbar

class DataPlotPanel(wx.Panel):
    """Class for plotting data from FitResultTable. Implements
    'observer' for FitResultTable."""
    
    colors = 'kbgcmr'
    """@ivar: list of colors for multiple plots"""

    colors_rgba = matplotlib.colors.ColorConverter().to_rgba_array(colors)

    markers = 'os^vd'
    """@ivar: list of marker styles for plots with multiple data"""

    linestyles = ['-', '--', '-.', ':']

    ID_AutoscaleX = wx.NewId()
    ID_AutoscaleY = wx.NewId()
    ID_Save       = wx.NewId()

    def __init__(self, parent, name = ''):
        """create panel and axes"""
        self.name = name
        wx.Panel.__init__(self, parent, size = (400, 300))
        self._create_canvas(parent)
        self.axs = []

    def _create_canvas(self, parent):
        self.fig = pylab.Figure()
        self.canvas = FigureCanvas(self, -1, self.fig)
        #self.canvas = FigureCanvas(parent, -1, self.fig)

        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.SetSizer(self.sizer)
        self._create_toolbar()
        tw, th = self.toolbar.GetSizeTuple()
        fw, fh = self.canvas.GetSizeTuple()
        self.toolbar.SetSize(wx.Size(fw, th))
        self.toolbar.update() #??
        self.sizer.Add(self.toolbar, 0, wx.LEFT | wx.EXPAND)

        self.Fit()
        
        self.canvas.mpl_connect('pick_event', self.onpick)

    def _create_toolbar(self):
        """create standard toolbar and add tools for toggling
        automatic scaling"""
        
        self.toolbar = NavigationToolbar2Wx(self.canvas)

        self.toolbar.DeleteTool(self.toolbar._NTB2_SAVE)
        self.toolbar.DeleteTool(self.toolbar._NTB2_SUBPLOT)
        

        self.toolbar.AddLabelTool(self.ID_Save,
                                  'Save',
                                  wx.Bitmap(os.path.join(settings.bitmappath,
                                                         'save.png'),
                                            wx.BITMAP_TYPE_PNG),
                                  shortHelp = "Save Image",
                                  longHelp = "Save Image to File")

        self.toolbar.AddSeparator()
        self.toolbar.AddCheckTool(self.ID_AutoscaleX,
                                  wx.Bitmap(os.path.join(settings.bitmappath,
                                                         'Xscale.png'),
                                            wx.BITMAP_TYPE_PNG),
                                  shortHelp = 'Autoscale X',
                                  longHelp = 'automatic scaling of X-axis')
        self.toolbar.AddCheckTool(self.ID_AutoscaleY,
                                  wx.Bitmap(os.path.join(settings.bitmappath,
                                                         'Yscale.png'),
                                            wx.BITMAP_TYPE_PNG),
                                  shortHelp = 'Autoscale Y',
                                  longHelp = 'automatic scaling of Y-axis')

        wx.EVT_TOOL(self, self.ID_Save, self.OnSave)
        wx.EVT_TOOL(self, self.ID_AutoscaleX, self.OnAutoscaleX)
        wx.EVT_TOOL(self, self.ID_AutoscaleY, self.OnAutoscaleY)

        self.toolbar.ToggleTool(self.ID_AutoscaleX, True)
        self.toolbar.ToggleTool(self.ID_AutoscaleY, True)

        self.autoscaleX = True
        self.autoscaleY = True

        self.toolbar.Realize()

    def draw(self):
        self.fig.canvas.draw()
        
    def update(self, subject):

        """update status by requesting and storing data from calling
        subject (given as argument), using the plotdata property"""

        #store pointer to subject (table) which sends message that
        #data has changed, needed for displaying mouse selection of
        #data
        self.subject = subject

        #get data from subject
        xdata, ydatas, d = subject.plotdata


        if d['xcol']: # and len(d['ycols']) > 0: #note xcol could be None
            self.datax = xdata
            self.datasy = ydatas
            self.dataid = d['yid']
            self.datagidx = d['gidx']
            self.datamask = d['masked']
            self.datadict = d

            self.refresh()

        else:
            #TODO: something wrong here!
            if False:
                #no valid data: escape
                self.datax = []
                self.datasy = []
                self.dataid = []
                self.datagidx = []
                self.datamask = []
                #TODO: fix for group data
                return


    def prepare_axes(self, n):
        """create n subplots (axes) by adding or removing subplots
        (axes), try to keep existing limits"""
        
        nold = len(self.axs)

        if n == nold:
            return

        #save old limits
        lims = [(a.get_xlim(), a.get_ylim()) for a in self.axs]

        #delete old axes
        for a in self.axs: self.fig.delaxes(a)
        self.axs = []

        #create new axes
        for k in range(n):
            if k>0:
                ax = self.fig.add_subplot(n,1,k+1, sharex = self.axs[0])
            else:
                ax = self.fig.add_subplot(n,1,k+1)
                
            if k>nold-1:
                #newly appended axis

                #need x autoscaling if no old axis exists
                if nold == 0:
                    ax.cam_need_autoscalex = True
                else:
                    ax.set_xlim(lims[0][0])
                    ax.cam_need_autoscalex = False
                
                #appended axes need autoscaling for first plot                    
                ax.cam_need_autoscaley = True
                
            else:
                ax.set_xlim(lims[k][0])
                ax.cam_need_autoscalex = False
                ax.set_ylim(lims[k][1])
                ax.cam_need_autoscaley = False
                
            ax.set_picker(True)
            ax.grid(True)
            self.axs.append(ax)


    def refresh(self):
        """refresh plot based on saved status"""

        if not self.canvas.IsShownOnScreen():
            return

        #plot data, store handles
        self.phs = []

        self.prepare_axes( len(self.datasy) )


        #do empty plots for group values legend
        #save axes limits
        xlimits = self.axs[0].get_xlim()
        ylimits = self.axs[0].get_ylim()
        self.axs[0].set_autoscale_on(False)
        gvals = self.datadict['gvals']
        sh = []
        st = []
        for k, gval in enumerate(gvals):
            h = self.axs[0].plot([], [], color = 'k', 
                                 linestyle = self.linestyles[k%len(self.linestyles)],
                                 )
            sh.append(h[0])
            st.append('%g'%gval)

        #restore limits
        self.axs[0].set_xlim(xlimits)
        self.axs[0].set_ylim(ylimits)

        # prepare 2 row of fields for saving averaged data of 2 subplots
        self.avgdatafields = []
        
        self.avg_plotted_data = []
        for k in range(len(self.datasy)):
            h = self.do_plot(self.axs[k],
                             self.datax, self.datasy[k],
                             self.datadict['ycollabels'][k],
                             self.datadict['ycols'][k],
                             groupidx = self.datagidx,
                             dataid = self.dataid,
                             showtitle = (k==0),
                             showxlabel = (k+1 == len(self.axs)),
                            )
            self.phs.append(h)
            
            fieldsrow = [self.datadict.get('xcollabel')]
            ylabels = self.datadict['ycollabels'][k]
            for icolumn in range(len(ylabels)):
                yfield = ylabels[icolumn]
                groups = numpy.unique(self.datagidx)
                if len(groups)>1:
                    yfield = ""
                    for gidx in groups:
                        if gidx<0:
                            pass
                        else:
                            fieldsrow.append(ylabels[icolumn]+"_"+st[gidx])
                            fieldsrow.append(ylabels[icolumn]+"_"+st[gidx]+"_err")
                else:
                    fieldsrow.append(ylabels[icolumn])
                    fieldsrow.append(ylabels[icolumn]+"_err")
            
            self.avgdatafields.append(fieldsrow)

        if sh:
            import matplotlib.legend
            l = matplotlib.legend.Legend(self.fig, sh, st, loc = 'upper left')
            self.axs[0].add_artist(l)
            #self.axs[0].legend(sh, st, loc = 'best')

        self.do_autoscale()
        
        #self.draw()
        #self.axs[0].get_xaxis().set_picker(True)

    def do_plot(self, ax, datax, datay, ylabels, ycols, groupidx, dataid, showtitle = True, showxlabel = True):
        """do plot in single axis. do not change axes limits"""

        if len(datax) == 0 or len(datay) == 0:
            print "no data, no plot"
            ax.clear()
            return #TODO: ????

        #create a masked array to store averaged data
        unique_datax = numpy.ma.array(numpy.unique(datax))
        avg_p_data = numpy.ma.vstack((unique_datax,numpy.ma.masked_all((2,len(unique_datax)))))
        avg_p_data = avg_p_data.transpose()
        
        #save axes limits
        xlimits = ax.get_xlim()
        ylimits = ax.get_ylim()

        ##workaround for 0.98.3, xlimits return view of internal data,
        ##need to copy data
        #if isinstance(xlimits, numpy.ndarray):
        #    xlimits = xlimits.copy()

        #if isinstance(ylimits, numpy.ndarray):
        #    ylimits = ylimits.copy()

        #clear and initialize axis
        ax.clear()
        ax.hold(True)
        ax.set_autoscale_on(False)

        handles = []

        for col in range(datay.shape[1]):
            try:
                groups = numpy.unique(groupidx)
                for gidx in groups:
                    if gidx<0 and len(groups)>1:
                        pass
                    else:
                        sel = groupidx == gidx
                        dysel = datay[sel, col]
                        #
                        dysel.shape = (-1, 1)
                        xmed, ymed, yerr = avgdata(datax[sel], dysel)
                        valid = ~ymed.mask.ravel()
                        #
                        ax.plot(xmed[valid],
                                ymed[valid],
                                color = self.colors[ycols[col]%len(self.colors)],
                                linestyle = self.linestyles[gidx%len(self.linestyles)],
                                )
                        ### fill the masked array of averaged data
                        for idatax in range(len(unique_datax)):
                            for ixmed in range(len(xmed[valid])):
                                if xmed[valid][ixmed] == avg_p_data[idatax,0]:
                                    avg_p_data[idatax,-2] = ymed[valid][ixmed]
                                    avg_p_data[idatax,-1] = yerr[valid][ixmed]
                        avg_p_data = numpy.ma.hstack((avg_p_data,numpy.ma.masked_all((len(unique_datax),2))))
                        
                        
                if False: #True:
                    #plot all data
                    h = ax.plot(datax, datay[:,col], 'o',
                               label = ylabels[col],
                               marker = self.markers[ycols[col]%len(self.markers)],
                               color = self.colors[ycols[col]%len(self.colors)],
                               picker = 3.0,
                               )
                    handles.extend(h)

                    #plot masked data
                    #TODO: error if all values of y are not valid
                    xd = datax[self.datamask]
                    yd = datay[self.datamask, col]

                    try:
                        ax.plot(datax[self.datamask],
                               datay[self.datamask, col],
                               'o',
                               label = '_nolegend_',
                               #marker = self.markers[ycols[col]%len(self.markers)],
                               markersize = 9,
                               markeredgecolor = 'b',
                               markerfacecolor = 'w',
                               alpha = 0.5,
                               )
                    except StandardError, e:
                       print "Error plotting masked data"
                       pass
                    
                else:
                    #use scatter to plot masked and non masked data
                    #points, using scatter. Note: problems with alpha
                    #handling in scatter need fixes, see
                    #sandbox/testscatter
                    cc = self.colors_rgba[[ycols[col]%len(self.colors)]*len(datax)]
                    cc[:,3] = 1.0 - self.datamask*0.7
                    #cc[dataid[-1],3] = 0.5

                    #TODO: emphasize last data point

                    h = ax.scatter(datax, datay[:,col],
                                   s = 15*groupidx+30 + 15,
                                   c = cc,
                                   marker = self.markers[ycols[col]%len(self.markers)],
                                   picker = 3.0,
                                   label = ylabels[col],
                                )
                    handles.extend([h])
                    
            except StandardError, e:
                print "Error in plotting data:", e
                raise
                return []
        
        # append to list of masked arrays containing avg data shown in plot
        # one entry of the list for each subplot
        self.avg_plotted_data.append(avg_p_data)
        
        #set axis labels, title, ...
        if showtitle:
            ax.set_title(label = self.datadict.get('name', ''), picker = True)
        else:
            ax.set_title(label = '')

        if showxlabel:
            ax.set_xlabel(self.datadict.get('xcollabel'), picker = 3.0)
        else:
            ax.set_xlabel('')

        #legend
        ax.legend(numpoints = 1, loc = 'best')

        #restore limits
        ax.set_xlim(xlimits)
        ax.set_ylim(ylimits)

        return handles
    
    def onpick(self, event):
        """handle pick events. If data point is hit, select it in the
        table (shift appends to current selection)"""

        if isinstance(event.artist, pylab.Line2D) or \
               isinstance(event.artist, matplotlib.collections.RegularPolyCollection):

            for subplotnr, handles in enumerate(self.phs):
                if event.artist in handles:
                    col = handles.index(event.artist)
                    break
            else:
                plotnr = None
                return
            
            line = event.artist
            ind = event.ind
            mm = (~self.datasy[subplotnr].mask[:, col]).nonzero()[0][ind]

            row = self.dataid[mm]
            self.subject.GetView().SelectRow(row, 
                                            (event.mouseevent.key == 'shift'))
            self.subject.GetView().MakeCellVisible(row, 0)
 
        elif isinstance(event.artist, matplotlib.axis.XAxis):
            print "you were hitting the xaxis"

        elif event.artist is self.axs[-1].get_xaxis().get_label():
            print "you were hitting the x label"
            print event.artist


    def OnAutoscaleX(self, event):
        self.autoscaleX = event.IsChecked()
        if self.autoscaleX:
            self.do_autoscale()

    def OnAutoscaleY(self, event):
        self.autoscaleY = event.IsChecked()
        if self.autoscaleY:
            self.do_autoscale()

    def do_autoscale(self):
        for ax in self.axs:
            ax.set_autoscale_on(True)
            ax.autoscale_view(scalex = self.autoscaleX or ax.cam_need_autoscalex,
                              scaley = self.autoscaleY or ax.cam_need_autoscaley)
            ax.cam_need_autoscaley = False
            ax.cam_need_autoscalex = False

        self.fig.canvas.draw()


    def OnSave(self, event):

        # Fetch the required filename and file type.
        filetypes, exts, filter_index = self.canvas._get_imagesave_wildcards()

        default_file = time.strftime("%Y%m%d") + '-' + self.name # + self.canvas.get_default_filetype()

        dlg = wx.FileDialog(self.canvas,
                            #_parent,
                            "Save to file", "",
                            default_file,
                            filetypes,
                            wx.SAVE|wx.OVERWRITE_PROMPT|wx.CHANGE_DIR)
        
        dlg.SetFilterIndex(filter_index)
        
        #
        #
        if dlg.ShowModal() == wx.ID_OK:
            dirname  = dlg.GetDirectory()
            filename = dlg.GetFilename()
            format = exts[dlg.GetFilterIndex()]
                        
            #Explicitly pass in the selected filetype to override the
            # actual extension if necessary
            try:
                self.canvas.print_figure(
                    str(os.path.join(dirname, filename)), format=format)
                
            except Exception, e:
                toolbar.error_msg_wx(str(e))
            
            self.save_avgdata_csv(dirname, filename)                
    
    def save_avgdata_csv(self, dirname, filename):
        """save data in comma separated format."""

        for isubplot in range(len(self.avg_plotted_data)):
            csvfilename = str(os.path.join(dirname, filename))
            # remove the extension
            csvfilename = csvfilename.rsplit(".")[0]
            csvfilename = csvfilename+'_'+str(isubplot)+'.csv'
        
            #don't save last two columns
            data   = self.avg_plotted_data[isubplot][:,:-2]

            with open(csvfilename, 'wb') as f:
                import csv
                writer = csv.writer(f)
                writer.writerow(self.avgdatafields[isubplot])
                for row in data:
                    r = [entry.encode('latin_1') if type(entry) is types.UnicodeType else entry for entry in row]
                    writer.writerow(r)
    
    
class DataTableApp(wx.App):

    def OnInit(self):
        self.frame = wx.Frame(None, title = "Data plot", size = (800, 700))
        self.plotpanel = DataPlotPanel(self.frame, 'test')
        self.gridpanel = FitResultTableGrid.GridPanel(self.frame)

        self.table = self.gridpanel.grid.Table
        self.grid  = self.gridpanel.grid

        self.statusbar = StatusBarWx(self.frame)
        self.plotpanel.toolbar.set_status_bar(self.statusbar)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.plotpanel, 1, wx.EXPAND)
        sizer.Add(self.gridpanel, 1, wx.EXPAND)
        sizer.Add(self.statusbar, 0, wx.LEFT | wx.TOP | wx.EXPAND)
        self.frame.Sizer = sizer
        self.frame.Show(True)

        self.create_data()

        #self.plotpanel.set_table(self.table)
        self.table.add_observer(self.plotpanel)

        #self.table.
        #self.plotpanel.update()

        #self.Bind(wx.EVT_ERASE_BACKGROUND, self.OnEraseBackground)
        return True

    def OnEraseBackground(self, event):
        print "Erase App!"
        event.Skip()


    def create_data(self):
        self.table.AppendRows(6)
        
        self.table.SetValueNamed(0, 'N K', 100)
        self.table.SetValueNamed(0, 'N Na',   10)
        self.table.SetValueNamed(0, 'user',   1)
        self.table.SetValueNamed(0, 'user2',  1)

        self.table.SetValueNamed(1, 'N K', 150)
        self.table.SetValueNamed(1, 'N Na',   20)
        self.table.SetValueNamed(1, 'user',   2)
        self.table.SetValueNamed(1, 'user2',  2)

        self.table.SetValueNamed(2, 'N K', 170)
        self.table.SetValueNamed(2, 'user',   3)
        self.table.SetValueNamed(2, 'user2',  1)
        
        self.table.SetValueNamed(3, 'N K', 200)
        self.table.SetValueNamed(3, 'N Na',   30)
        self.table.SetValueNamed(3, 'user2',  2)

        self.table.SetValueNamed(4, 'N K', 250)
        self.table.SetValueNamed(4, 'N Na',   40)
        self.table.SetValueNamed(4, 'user',   5)
        self.table.SetValueNamed(4, 'user2',  1)

        for row, value in enumerate([10, 20, 30, 30, 50]):
            self.table.SetValueNamed(row, 'sx Na', value)
            self.table.SetValueNamed(row, 'sxerr Na', value**0.5)

        #self.grid.SetColumnSelection([0, 2, 13, 24, 28, 29, 31, 32])
        self.grid.SetColumnSelection([self.table.colname_to_raw('FileID'),
                                      self.table.colname_to_raw('N Na'),
                                      self.table.colname_to_raw('N K'),
                                      self.table.colname_to_raw('sx Na'),
                                      self.table.colname_to_raw('sxerr Na'),
                                      self.table.colname_to_raw('dynamic'),
                                      self.table.colname_to_raw('user'),
                                      self.table.colname_to_raw('user2'),
                                      self.table.colname_to_raw('Omit'),
                                      self.table.colname_to_raw('Remark'),
                                      ])

        self.table.SetColMark(6, 'X')
        self.table.SetColMark(1, 'Y1')
        self.table.SetColMark(2, 'Y1')
        self.table.SetColMark(2, 'Y2')

        self.table.dynamic_expressions[0] = 'N_Na + N_K'


def test():
    #gui = DataApp()
    gui = DataTableApp(redirect = False)
    gui.MainLoop()
    return gui


if __name__ == '__main__':
    gui = test()
    
