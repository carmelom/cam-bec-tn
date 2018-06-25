#!/usr/bin/python
#-*- coding: latin-1 -*-

import os.path
import wx
import numpy as np
from scipy.optimize import curve_fit

import csv

from scipy.integrate import trapz, quad
from scipy.signal import fftconvolve

def convolution(f, g, x, sigma, relerr = 1e-5):
    """
    Performs the numerical convolution of f by g, evaluted at the datapoints x. It is assumed that at least one function presents an exponential decay, on caracteristic lenght sigma.
    
    Parameters
    ----------
    f, g : function
        Functions taking as parameters and returning a 1d-ndarrays
    x : 1d-ndarray
        Data points where the convolution is evaluated
    sigma : float
        Caracteristic width of the thinner distribution
    """
    xLim = - sigma*np.log(relerr)
    dx = xLim * relerr
    xp = np.arange(x.min() - xLim, x.max() + xLim, dx)
    fp = dx * fftconvolve(f(xp), g(xp), mode='same')
    return np.interp(x, xp, fp)

# convolution(lambda x: 1/((1+(x/fr)**2))*np.sin(17e-3*np.pi*np.sqrt(fr**2+x**2))**2, lambda x: 1/(np.sqrt(2*np.pi)*fd)*np.exp(-.5*(x/fd)**2), x-fc, 30)
# fc, fr, fd
# 550, 22, 10

class FitManager():
    
    NbPtsFit = 256
    nameColSave = ['name', 'nameX', 'nameY', 'nameG', 'nameGVal', 'enableG', 'fct', 'params', 'fitParams', 'relErr', 'errorsParams', 'omit',]
    
    def __init__(self, dataPanel):
        
        self.dataPanel = dataPanel
        
        self.ID_fitPanel            = wx.NewId()
        self.ID_RightClickReloadY1  = wx.NewId()
        self.ID_RightClickReloadY2  = wx.NewId()
        self.ID_RightClickDelete    = wx.NewId()
        
        self.fitPanel = FitPanel(dataPanel, self)
    
        self.fits = {}
        
        self.compatibleFits = {} # dict, bool
        self.XYGParams = False
        
        self.previewActive = False
        
    def _getDataY(self, entries):        
        YChoiceName = self.XYGChoices['Y'][entries['Y']]
        for i,yLabels in enumerate(self.dataPanel.datadict['ycollabels']):
            for j,yLabel in enumerate(yLabels):
                if YChoiceName == yLabel:
                    ydata = self.dataPanel.datasy[i][:,j]
                    yid = self.dataPanel.ycol[i][j]
        return ydata, yid
    
    def _getEntries(self):
        entries = {}
        for name, field in self.dataPanel.fields.items():
            if isinstance(field, wx.TextCtrl):
                entries[name] = field.GetLineText(0)
            elif isinstance(field, wx.Choice):
                entries[name] = field.GetSelection()
            elif isinstance(field, wx.CheckBox):
                entries[name] = field.GetValue()
        return entries
    
    
    def doFit(self, evt):
        if self.XYGParams:
            
            entries = self._getEntries()
            
            compatible = True
            enableG = False
                    
            # For 'G' and 'GVal', the 1st elt (num 0) is '/' (no entry)
            if bool(entries['G']) !=  bool(entries['GVal']):               # Or G and GVal have no entry, or they have
                    compatible = False
            elif entries['G']:
                enableG = True
                    
            if len(entries['name']) == 0:
                compatible = False
                
            if compatible:
                      
                xdata = self.dataPanel.datax
                ydata, yid = self._getDataY(entries)
                
                mask = ~(np.isnan(xdata) | np.isnan(ydata))
                
                if enableG:
                    mask &= (self.dataPanel.datagidx + 1 == entries['GVal'])
                 
                if entries['relErr']:
                    yerr = np.abs(ydata[mask])
                else:
                    yerr = None
                    
                try:                    
                    fct = eval("lambda x,{params} : {fct}".format(params=entries['params'], fct=entries['fct']))
                    popt, mcov = curve_fit(fct, xdata[mask], ydata[mask], sigma=yerr, p0=eval(entries['initParams']))
                except RuntimeError:
                    print "Error: The fit doesn't converge"
                    return
                except Exception as e:
                    print "Error:", e
                    return
                    
                nameFit = entries['name']
                fit = self.fits[nameFit] = entries.copy()
                
                fit['mask'] = mask
                fit['fctNum'] = fct
                fit['fitParamsNum'] = popt
                fit['fitParams'] = ''.join(["%.3e, "%p for p in popt])
                fit['errorsParamsNum'] = np.diag(mcov)
                fit['errorsParams'] = ''.join(["%.2e, "%e for e in self.fits[nameFit]['errorsParamsNum']])
                fit['xdata'] = xdata
                
                for name, choices in self.XYGChoices.items():
                    fit["name"+name] = choices[entries[name]]
                
                fit['xid'] = self.dataPanel.xcol[0]
                fit['yid'] = yid
                fit['gid'] = self.dataPanel.gcol[0] if enableG else None
                fit['relErr'] = entries['relErr']
                fit['relErrStr'] = "1" if entries['relErr'] else "0"
                fit['enableG'] = enableG
                fit['omit'] = ''   # bool cells: '' if False, '1' if True
                
                self.compatibleFits[nameFit] = True
                self.fitPanel.update(self.fits)
                self.dataPanel.refresh()
    
    def refresh(self, XYGChoices):
        self.XYGChoices = XYGChoices
        self.XYGParams = True
        for choices in self.XYGChoices.values():
            if not choices:
                self.XYGParams = False
        if self.XYGParams:
            for fitName in self.fits.keys():
                compatible = True
                for name, choices in self.XYGChoices.items():
                    if self.fits[fitName]["name"+name] not in choices:
                        compatible = False
                self.compatibleFits[fitName] = compatible
                
            self.fitPanel.update(self.fits)
            self.firstPreview = True
    
    def _SetColMark(self, table, colId, mark):
        return table.SetColMark(table.colsel.index(colId), mark)
    
    def reload(self, nameFit, markY):
        fit = self.fits[nameFit]
        fitResultdataTableGrid = self.dataPanel.gridTable.GetFitResultDataTableGrid()
        table = fitResultdataTableGrid.Table
        
        self._SetColMark(table, fit['xid'], 'X')
        self._SetColMark(table, fit['yid'], markY)
        if fit['enableG']:
            self._SetColMark(table, fit['gid'], 'G')
            
        fitResultdataTableGrid.Refresh()
        
        for name, field in self.dataPanel.fields.items():
            if isinstance(field, wx.TextCtrl):
                field.ChangeValue(fit[name])
            elif isinstance(field, wx.Choice):
                field.SetSelection(fit[name])
            elif isinstance(field, wx.CheckBox):
                field.SetValue(fit[name])
    
    def delete(self, nameFit):
        del self.fits[nameFit]
        del self.compatibleFits[nameFit]
        self.fitPanel.update(self.fits)
        self.dataPanel.refresh()
        
    def omit(self, nameFit, state):
        self.fits[nameFit]['omit'] = state
        self.dataPanel.refresh()
        
    def _getChoices(self):
        self.XChoices = self.dataPanel.XChoices
        self.YChoices = self.dataPanel.YChoices
        self.GChoices = self.dataPanel.GChoices
    
    # Each time a fit is performed or the axis modified, all the fits are deleted and reploted
    def plot(self):
        self.previewActive = False
        nbColors = len(self.dataPanel.colors)
        for ax, yLabels, yCol in zip(self.dataPanel.axs, self.dataPanel.datadict['ycollabels'], self.dataPanel.ycol):
            
            for nameFit, fit in self.fits.items():
                if self.compatibleFits[nameFit] and not fit['omit'] and fit['nameY'] in yLabels:
                    
                    X = np.linspace(fit['xdata'].min(), fit['xdata'].max(), self.NbPtsFit)
                    Y = fit['fctNum'](X, *fit['fitParamsNum'])
                    if not isinstance(Y, np.ndarray):
                        Y *= np.ones(X.shape)
                    
                    yNum = yCol[yLabels.index(fit['nameY'])]
                    ax.plot(X, Y, marker=None, linestyle='-', color=self.dataPanel.colors[yNum%nbColors])
    
    # Enable the setting of init parameters by the preview of the plotted function
    def preview(self, evt):
        if self.XYGParams:
            xdata = self.dataPanel.datax
            entries = self._getEntries()
            
            xMin, xMax = self.dataPanel.xLim
            
            for i,(ax, yLabels) in enumerate(zip(self.dataPanel.axs, self.dataPanel.datadict['ycollabels'])):
                yMin, yMax = self.dataPanel.yLim[i]
                
                if self.XYGChoices['Y'][entries['Y']] in yLabels:
                    
                    if self.previewActive:
                        del ax.lines[-1]
                        self.previewActive = False
                    
                    X = np.linspace(xMin, xMax, self.NbPtsFit)
                    try:
                        fct = eval("lambda x,{params} : {fct}".format(params=entries['params'], fct=entries['fct']))
                    
                        p0 = eval(entries['initParams'])
                        if not isinstance(p0, tuple):
                            p0 = (p0,)
                            
                        Y = fct(X, *p0)
                        if not isinstance(Y, np.ndarray):
                            Y *= np.ones(X.shape)
                        yMin = min(yMin,Y.min())
                        yMax = max(yMax,Y.max())
                        
                        ax.plot(X, Y, marker=None, linestyle=(0,(5,10)), linewidth=3, color='silver')
                        self.previewActive = True
                    except:
                        print "Error: invalid parameters for preview" 
                  
                ax.set_xlim(xMin, xMax)
                ax.set_ylim(yMin, yMax)
                
            self.dataPanel.draw()
    
    def save(self, dirname, filename):
        
        path = str(os.path.join(dirname, filename)).rsplit(".")[0] + "-fit.csv"
        
        with open(path, 'wb') as f:
            writer = csv.writer(f, delimiter=';')
            writer.writerow(self.nameColSave)
            
            for fit in self.fits.values():
                saveFitsRow = []
                for name in self.nameColSave:
                    saveFitsRow.append(fit[name])
                writer.writerow(saveFitsRow)          
                    
                #r = [entry.encode('latin_1') if type(entry) is types.UnicodeType else entry for entry in row]


        
              
class FitPanel(wx.grid.Grid):
    
    gray = wx.Colour(220, 220, 220)
    red = wx.Colour(230,10,60)
    entriesNames = 'name', 'nameX', 'nameY', 'nameG', 'nameGVal',      'fct',     'params',     'fitParams',  'relErrStr', 'errorsParams', 'omit'
    colNames =     'name',     'X',     'Y',     'G',     'GVal', 'function', 'parameters', 'fitted values',    'rel err',       'errors', 'omit'
    colWidths =       100,      50,      50,      50,         50,        190,          110,             250,           50,            250,     50
    iOmit = colNames.index('omit')
    
    def __init__(self, dataPanel, fitManager):
    
        """self.data contains the informations relating to the fits, in the following order:
           name, X data, Y data, G group by, fonction definition, parameters (..,..,..), and parameters initial values"""
                    
        self.fitManager = fitManager
        wx.grid.Grid.__init__(self, dataPanel, self.fitManager.ID_fitPanel, name="fit list")
        
        self.data = {}
        self.sortedNames = []
        self.rowRightclick = 0
        
        self.CreateGrid(0, len(self.colNames), selmode=self.SelectRows)
        
        for i, (colName, colWidth) in enumerate(zip(self.colNames, self.colWidths)):
            self.SetColLabelValue(i, colName)
            self.SetColSize(i, colWidth)
        
        self.SetColFormatCustom(self.iOmit, wx.grid.GRID_VALUE_BOOL)
        
        wx.grid.EVT_GRID_LABEL_RIGHT_CLICK(self, self._rightClickRow)
        wx.grid.EVT_GRID_CELL_RIGHT_CLICK(self, self._rightClickRow)
        wx.grid.EVT_GRID_CELL_CHANGING(self, self._clickOmit)
        
        self.Bind(wx.EVT_MENU, self._clickReloadY1, id=self.fitManager.ID_RightClickReloadY1)
        self.Bind(wx.EVT_MENU, self._clickReloadY2, id=self.fitManager.ID_RightClickReloadY2)
        self.Bind(wx.EVT_MENU, self._clickDelete, id=self.fitManager.ID_RightClickDelete)
    
    def update(self, fitResults):
        self.sortedNames = np.sort(fitResults.keys())
        
        self.ClearGrid()
        nbRows = self.GetNumberRows()
        if nbRows:
            self.DeleteRows(numRows=self.GetNumberRows())
        self.AppendRows(len(fitResults))
        for i,nameFit in enumerate(self.sortedNames):   
            for j,name in enumerate(self.entriesNames):
                self.SetCellValue(i,j,fitResults[nameFit][name])
                if not self.fitManager.compatibleFits[nameFit]:
                    self.SetCellBackgroundColour(i,j,self.gray)
                    self.SetCellTextColour(i,j,self.red)

    def _rightClickRow(self, evt):
        self.rowRightclick = evt.Row
        if self.rowRightclick >= 0: #right click on a row
            self.SelectRow(self.rowRightclick)
            
            menu = wx.Menu()
            menu.Append(self.fitManager.ID_RightClickReloadY1,
                              "Reload Y1")
            menu.Append(self.fitManager.ID_RightClickReloadY2,
                              "Reload Y2")
            menu.Append(self.fitManager.ID_RightClickDelete,
                              "Delete")
            self.PopupMenu(menu)
            menu.Destroy()
        
    def _clickReload(self, evt, markY):
        if self.rowRightclick >= 0:
            self.fitManager.reload(self.sortedNames[self.rowRightclick], markY)
            
    def _clickReloadY1(self, evt):
        self._clickReload(evt, 'Y1')
            
    def _clickReloadY2(self, evt):
        self._clickReload(evt, 'Y2')
        
    def _clickDelete(self, evt):
        if self.rowRightclick >= 0:
            self.fitManager.delete(self.sortedNames[self.rowRightclick])
    
    def _clickOmit(self, evt):
        name = self.sortedNames[evt.Row]
        state = evt.GetString()
        self.fitManager.omit(name, state)
