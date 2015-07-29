'''
*********************************************************************
gkPluginInterface.py

Defining the user interactions of the plugin's Graphical User Interface (GUI).
GUI component classes are defined in form.py
*********************************************************************
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

   Begin     : 2012-07-01
   Copyright : (C) 2012 Go Sato 
   Email     : go.sato@zaq1.net
*********************************************************************
'''

from PyQt4.QtCore import *
from PyQt4.QtGui import *
from qgis.core import *
from form import Ui_Dialog
from gausskruger import transform2GaussKruger
from wkt import crackWKT
import resources_rc
import os
import re


class gkPluginInterface:

  def __init__(self, iface):
      # Save reference to the QGIS interface
      self.iface = iface
      self.canvas = iface.mapCanvas()
      self.curDir = QString(os.getenv('USERPROFILE') + "\\My Documents")
      self.fileDiag = QFileDialog()
    
  
  def debugMsg(self, msg):
    QMessageBox.information(self.iface.mainWindow(), "Message", QString(str(msg)))

  def initGui(self):
      # Create action that will start plugin
      self.action = QAction(QIcon(":/transmerc_icon2.png"), "TransMerc", self.iface.mainWindow())
      
      # connect the action to the run method
      QObject.connect(self.action, SIGNAL("triggered()"), self.run)

      # Add toolbar button and menu item
      self.iface.addPluginToMenu("&TransMerc", self.action)
      

  def unload(self):
      # Remove the plugin menu item and icon
      self.iface.removePluginMenu("&TransMerc",self.action)

  # run
  def run(self):
    self.win = QDialog()
    self.win.ui = Ui_Dialog()
    self.win.ui.setupUi(self.win)
    QObject.connect(self.win.ui.transformBtn, SIGNAL("clicked()"), self.close)
    QObject.connect(self.win.ui.selectQGisLayer, SIGNAL("currentIndexChanged(int)"), self.displayLayerProjection)
    QObject.connect(self.win.ui.fileDialogBtn, SIGNAL("clicked()"), self.fileSelect)
    QObject.connect(self.win.ui.saveFilePath, SIGNAL("textEdited(QString)"), self.enableTransformBtn)
    
    
    self.layerSelected = False
    self.winFilePathREGEX = re.compile(r'[a-zA-z]:/([^\\\*\?\"<>/]+/)*[^\\\*\?\"<>/]+')
    self.setupLayerList()
    self.displayLayerProjection()
    self.win.show()
    self.win.exec_()
    return

  def close(self):
    layer = self.getSelectedLayer()
    oLon = self.win.ui.originLongitudeField.value()
    oScl = self.win.ui.originScale.value()
    pBar = self.win.ui.progressBar
    falseN = self.win.ui.falseNorthingField.value()
    falseE = self.win.ui.falseEastingField.value()
    savePath = str(self.win.ui.saveFilePath.text())
    loadFile = self.win.ui.checkBox.isChecked()
    gk = transform2GaussKruger(layer, oLon, oScl, falseE, falseN, savePath, loadFile, pBar)
    gk.transformVLayerToGK()
    #self.debugMsg(layer.name())

  def setupLayerList(self):
    layers = QgsMapLayerRegistry.instance().mapLayers()
    for (key, layer) in layers.iteritems():
      if layer.isValid() and layer.type() == 0:
        flag = True
        if layer.dataProvider().crs().epsg() == 0:
          bbx = layer.extent()
          if bbx.xMaximum() > 180.0 or bbx.xMinimum() < -180.0 or \
             bbx.yMaximum() >  90.0 or bbx.yMinimum() <  -90.0:
            flag = False
        if flag:    
          self.win.ui.selectQGisLayer.addItem(layer.name(), QVariant(layer.getLayerID()))
  
  def getSelectedLayer(self):
    if self.win.ui.selectQGisLayer.count() > 0:
      i = self.win.ui.selectQGisLayer.currentIndex()
      layerID = self.win.ui.selectQGisLayer.itemData(i).toString()
      return QgsMapLayerRegistry.instance().mapLayer(layerID)
    else:
      return None
  
  def displayLayerProjection(self):
    layer = self.getSelectedLayer()
    if layer != None:
      crs = layer.dataProvider().crs()
      if crs.epsg() > 0:
        html = crackWKT(str(crs.toWkt())).exportDescription()
      else:
        html = '''
          <span style="color: red">Coordinate reference system (CRS) is not defined for this layer<br><br>
          Are you sure that this layer is WGS84?</span>
        '''
      self.layerSelected = True
    else:
      html = '<span style="color: red">There are no valid vector layers opened in this project!</span>'
    self.win.ui.layerProjectionDisplay.clear()
    qc = QTextCursor()
    self.win.ui.layerProjectionDisplay.appendHtml(QString(html))
    self.win.ui.layerProjectionDisplay.moveCursor(qc.Start)
    self.win.ui.layerProjectionDisplay.ensureCursorVisible()
  
  def fileSelect(self):
    filters = "ESRI Shapefile (*.shp *.SHP);;" + \
              "AutoCAD DXF (*.dxf *.DXF);;"    + \
              "Mapinfo (*.mif *.MIF);;"        + \
              "Geographic Markup Language (*.gml *.GML)"    
    path = self.fileDiag.getSaveFileName(None, "Save to", self.curDir, filters)
    if path != '':
      self.win.ui.saveFilePath.setText(path)
      r = str(path).split('/')
      r.pop()
      self.curDir = '/'.join(r)
    self.enableTransformBtn()
  
  def modifyCurDir(self, directory):
    self.curDir = directory
    
  
  def enableTransformBtn(self):
    btn = self.win.ui.transformBtn
    enabled = btn.isEnabled()
    selected = self.layerSelected
    saveFile = False
    path = self.win.ui.saveFilePath.text()
    m = self.winFilePathREGEX.match(path)
    if m and m.group(0) == path:
      fileName = (path.split('/'))[-1]
      fNamElems = fileName.split('.')
      extension = str(fNamElems[-1])
      if extension.lower() != 'shp' and extension.lower() != 'dxf' and \
         extension.lower() != 'mif' and extension.lower() != 'gml':
        pass
      else:
        for i in range(len(fNamElems)-1):
          if len(fNamElems[i]) > 0:
            saveFile = True;  break
    
    if not enabled and selected and saveFile:
      btn.setEnabled(True)
    elif enabled and (not selected or not saveFile):
      btn.setEnabled(False)
        

if __name__ == "__main__":
    pass

