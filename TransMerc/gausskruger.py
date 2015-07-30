'''
*********************************************************************
gausskruger.py

Main functions for transforming the projection of QGIS layers to
Transverse Mercator using Gauss Kruger equations.
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

from qgis.core import *
from karney_gsl import TransverseMercatorExact
from pointBox import pointBox
from PyQt4.QtCore import QVariant
from PyQt4.QtGui import QApplication, QDialog, QMessageBox
import sys
import math

dummy = QgsGeometry()
Equator = dummy.fromWkt('LINESTRING(-181.0 0.0, 181.0 0.0)')
northHemi = dummy.fromWkt("POLYGON ((-181.0 0.0000000001, 181.0 0.0000000001, 181.0 91.0, -181.0 91.0, -181.0 0.0000000001))")
southHemi = dummy.fromWkt("POLYGON ((181.0 -0.0000000001, -181.0 -0.0000000001, -181.0 -91.0, 181.0 -91.0, 181.0 -0.0000000001))")


'''
Calculates the perpendicular distance from point (x, y) to a line
which pass through points (x1, y1) and (x2, y2)
Referenced by QGisPointToLineDistance()
References "math" library
'''
def pointToLineDistance(x, y, x1, y1, x2, y2):
  a = y1 - y2;  b = x2 - x1;  c = x1 * y2 - x2 * y1
  d = math.hypot(a, b)
  if d == 0:
    return math.hypot(x - x1, y - y1)
  else:
    return math.fabs(a*x + b*y + c) / d

'''
Calculates the perpendicular distance from point (qpt) to a line
which pass through points (qpt1) and (qpt2)
qpt, qpt1 and qpt2 are all QgsPoint type.
Referenced by transform2GaussKruger.transoformLineSection()
References qgis.core library and pointToLineDistance()
'''
def QGisPointToLineDistance(qpt, qpt1, qpt2):
  return pointToLineDistance(qpt.x(), qpt.y(), qpt1.x(), qpt1.y(), qpt2.x(), qpt2.y())

'''
transform2GaussKruger Class holds the parameters for Transverse Mercator
projection (TM) given by the users and provides means to transform the 
projections of input layers to TM.
Methods in this class references the following resources
  math         : Python standard library
  qgis.core    : Quantum GIS API
  karney_gsl.py: Gauss-Kruger calculation engine  
'''
class transform2GaussKruger():
  '''
  Constructor. the variables are as follows
    qgisLayer (QgsVectorLayer): 
      The layer selected by the user out of them loaded in QGIS
    originLongitude (double): 
      The longitude of the prime meridian given by the user
    originScale (double): 
      The scale factor
    falseEasting (double): 
      False easting given by the user
    falseNorthing (double): 
      False northing given by the user
    saveFilePath (string): 
      The file-path to save the output (transformed) map
    loadOutput (boolean): 
      Whether the output should be loaded into current project or not.
    QProgBar (QtGui.QProgressBar):
      The PyQt4 progress bar object to which the progress of transformation is shown    
  '''
  def __init__(self, qgisLayer, originLongitude, originScale, falseEasting, falseNorthing, saveFilePath, loadOutput, QProgBar):
    self.layer = qgisLayer
    self.QProgBar = QProgBar
    self.lon0 = originLongitude
    self.scale = originScale
    self.falseE = falseEasting
    self.falseN = falseNorthing
    self.saveFilePath = saveFilePath
    
    # Determine the file driver used for output from the extension in
    # the save-file path
    extension = (((self.saveFilePath.split('/'))[-1]).split('.'))[-1]
    if extension.lower() == 'dxf':
      self.fDriver = "DXF"
    elif extension.lower() == 'mif':
      self.fDriver = "MapInfo File"
    elif extension.lower() == 'gml':
      self.fDriver = "GML"
    else:
      self.fDriver = "ESRI Shapefile"
    
    self.WKBType = self.layer.wkbType() 
    self.loadAfterTransformation = loadOutput
    self.tolerance = 1000
    self.TrvMercator = TransverseMercatorExact(6378137.0, 6356752.314140, self.lon0, self.scale, False)
    self.bLon = math.degrees((1.0 - self.TrvMercator.e) * math.pi*0.5)
    self.provider = self.layer.dataProvider()
    self.epsg = self.provider.crs().authid()
    # QMessageBox.information(None, "DEBUG:", str(self.epsg)) 
    
    CRS1 = QgsCoordinateReferenceSystem()
    CRS1.createFromString(self.epsg)
    CRS2 = QgsCoordinateReferenceSystem()
    CRS2.createFromString('EPSG:4326')
    # CRS2 = QgsCoordinateReferenceSystem(4326)
    self.xform = QgsCoordinateTransform(CRS1, CRS2)
    
    ft = QgsFeature()
    # self.allAttribs = self.provider.attributeIndexes()
    # self.provider.select(self.allAttribs)
    feats = self.provider.getFeatures()
    self.vNum = 0
    while feats.nextFeature(ft):  # self.provider.nextFeature(ft):
      geom = ft.geometry()
      self.vNum += self.countVertices(geom)
    self.transformedPoints = 0
  
  '''
  Main routine for transforming the projection of Quantum GIS's Vector 
  Layer to Gauss-Kruger.
  '''
  def transformVLayerToGK(self):
    provider = self.provider  
    feat = QgsFeature()
    # self.provider.select(self.allAttribs)
    feats = self.provider.getFeatures()
  
    fields = self.provider.fields()
    writer = QgsVectorFileWriter(self.saveFilePath, "CP1250", fields, self.WKBType, None, self.fDriver)
  
    if writer.hasError() != QgsVectorFileWriter.NoError:
      print "Error when creating shapefile: ", writer.hasError()
  
    while feats.nextFeature(feat): # self.provider.nextFeature(feat):
      geom = feat.geometry()
      attr = feat.attributes()
      pb = pointBox(geom)
      self.transformGeomToWGS84(pb)
      
      splitFlag = False
      if pb.category != 'Point':
        splitFlag = self.checkEquatorCrossing(pb.toQGSGeometry())
        
      if splitFlag: 
        pb = self.splitGeometryByEquator(pb)
        
      self.latlon2XY(pb)  
      
      if splitFlag:
        pb = self.mergeGeometriesAtEquator(pb)
      self.applyScaleAndFalseNE(pb)  
      
      fet = QgsFeature()
      fet.setGeometry(pb.toQGSGeometry())

      fet.setAttributes(attr)
      writer.addFeature(fet)
      
    del writer
    if self.loadAfterTransformation:
      layer = QgsVectorLayer(self.saveFilePath, "layer_name_you_like", "ogr")
      QgsMapLayerRegistry.instance().addMapLayer(layer)


  '''
  Transforms all the coordinates of points in a pointBox object from
  (longitude, latitude) to (Easting, Northing) of Gauss-Kruger projection.
  All the points coordinates in pointbox must fall in the following domains
  X: [-180, 180] (longitude in decimal degrees)
  Y: [ -90,  90] (latitude in decimal degrees)
  
  Referenced by self.transformVLayerToGK()
  ''' 
  def latlon2XY(self, pointbox):
    polyBoundaries = []
    if pointbox.category != 'Point':
      for lv2 in pointbox.vector:
        for j in range(len(lv2)):
          lv1 = lv2[j];  newLv1 = [];
          for i in range(len(lv1) - 1):
            pt1 = lv1[i];  pt2 = lv1[i+1]
            segmentVertices = []
            null = None
            self.transoformLineSection(pt1, pt2, segmentVertices, null, null)
            if i < len(lv1) - 2:
              segmentVertices.pop()
            else:
              self.progressBar()
                  
            newLv1.extend(segmentVertices)
            self.progressBar()
            
          for i in range(len(newLv1)):
            p = newLv1[i]
            newLv1[i] = QgsPoint(p.x(), p.y())
           
          lv2[j] = newLv1
    else:
      for lv2 in pointbox.vector:
        for lv1 in lv2:
          for i in range(len(lv1)):
            pt = lv1[i]
            (x, y, g, k) = self.TrvMercator.Forward(pt.y(), pt.x())
            lv1[i] = QgsPoint(x, y)
            self.progressBar()
                    
    pointbox.checkWkbType()    

  '''
  Projects the line between lat-lon points qpt1 and qpt2 to
  X-Y plane using Gauss-Kruger projection. 
  This method approximates the curvature of the projected 
  line by dividing the input line at the midpoint and recursively 
  applying itself for each division.
  This method references the class parameter 'tolerance' as the 
  tolerance of approximation.
  
  Variables:
  qpt1 (QgsPoint):
    the start lat-lon point of the line segment
  qpt2 (QgsPoint):
    the end lat-lon point of the segment
  line (list):
    the list to store the constituent QgsPoint's of the output
    Users must give an empty list for this variable.
  qpt1x, qpt2x (QgsPoint):
    Only used for saving time for recursive calls.
    Users must give 'None' for these variables.
    
  Referenced by self.latlon2XY()      
  '''  
  def transoformLineSection(self, qpt1, qpt2, line, qpt1x, qpt2x):
    if qpt1x == None:
      (x, y, g, k) = self.TrvMercator.Forward(qpt1.y(), qpt1.x())    
      p1 = QgsPoint(x, y)    
    else:
      p1 = qpt1x
    
    if len(line)==0 or line[-1] != p1:
      line.append(p1)
    
    if qpt2x == None:
      (x, y, g, k) = self.TrvMercator.Forward(qpt2.y(), qpt2.x())    
      p2 = QgsPoint(x, y)
    else:
      p2 = qpt2x    
    
    midLon = (qpt1.x() + qpt2.x()) * 0.5
    midLat = (qpt1.y() + qpt2.y()) * 0.5
    (xm, ym, g, k) = self.TrvMercator.Forward(midLat, midLon)
    pm = QgsPoint(xm, ym)
  
    if QGisPointToLineDistance(pm, p1, p2) > self.tolerance:
      qptm = QgsPoint(midLon, midLat)
      self.transoformLineSection(qpt1, qptm, line, p1, pm)
      self.transoformLineSection(qptm, qpt2, line, pm, p2)
    else:
      line.append(p2)  
  
  
  '''
  Checks whether the given QgsGeometry object (polyGeom) would be 
  split into two at the equator by Gauss-Kruger projection.
  polyGeom must be a geometry with lat-lon coordinates.
  
  Referenced by self.transformVLayerToGK()   
  '''
  def checkEquatorCrossing(self, polyGeom):
    flag = False
    if polyGeom.crosses(Equator):
      bb = polyGeom.boundingBox()
      lonMin = self.normalisedLongitude(bb.xMinimum())
      lonMax = self.normalisedLongitude(bb.xMaximum())
      if lonMax > self.bLon or lonMin < -self.bLon or lonMax < lonMin:
        flag = True
    return flag  
  
  '''
  Splits the given pointBox object into two at the equator.
  To avoid confusions in the Gauss-Kruger transformation, the 
  cutting edge would be pushed away from the equator by 0.0000000001
  degrees.
  polyGeom must be defined with lat-lon coordinates.
  
  Referenced by self.transformVLayerToGK()  
  '''
  def splitGeometryByEquator(self, pointbox):
    if pointbox.category == 'Point':
      return pointbox
    else:
      geom = pointbox.toQGSGeometry()
      northPart = northHemi.intersection(geom)
      southPart = southHemi.intersection(geom)
      newPB = pointBox(northPart);
      southPB = pointBox(southPart);
      newPB.simplyUnite(southPB)
      return newPB
      
  '''
  First checks whether there are any polygons/polylines in the given 
  pointBox object which are sticking together at the equator.
  If yes, it tries to merge the adjacent geometries into a single
  geometry.
  This method typically works for polygons which surround the 
  singularity point of Gauss-Kruger equations (the point where the 
  equator starts to split toward north and south in Gauss-Kruger
  projection).
  
  Referenced by self.transformVLayerToGK()
  '''
  def mergeGeometriesAtEquator(self, pointbox):    
    if pointbox.category == 'Point':
      return pointbox    
    elif pointbox.category == 'Polyline':
      ln = len(pointbox.vector[0])
      if ln > 1:
        for i in range(ln - 1, 1, -1):
          pt1 = []
          pt1.append(pointbox.vector[0][i][ 0]);  
          pt1.append(pointbox.vector[0][i][-1]);
          connectThem = False
          if math.fabs(pt1[0].y()) < 0.0001 or math.fabs(pt1[1].y()) < 0.0001 :
            reverseI = None;  reverseJ = None
            for j in range(i-1, 0, -1):
              pt2 = []
              pt2.append(pointbox.vector[0][j][ 0])  
              pt2.append(pointbox.vector[0][j][-1])
              if math.fabs(pt2[0].y()) < 0.0001 or math.fabs(pt2[1].y()) < 0.0001 :
                for k in range(len(pt1)):
                  for l in range(len(pt2)):
                    if math.fabs(pt1[k].y()) < 0.0001 and math.fabs(pt2[l].y()) < 0.0001 and \
                       math.fabs(pt1[k].x() - pt2[l].x()) < 0.01:
                      connectThem = True
                      if k == 1:
                        pointbox.vector[0][i].reverse()
                      if l == 1:
                        pointbox.vector[0][j].reverse()  
                      break
                  if connectThem:
                    break
                  
                if connectThem:  
                  k = 0
                  for pt in pointbox.vector[0][i]:  
                    if k == 0:
                      x = (pt.x() + pointbox.vector[0][j][0].x()) * 0.5
                      pointbox.vector[0][j][0] = QgsPoint(x, 0.0)
                    else:  
                      pointbox.vector[0][j].insert(0, pt)
                    k += 1
                  break    
          if connectThem:
            pointbox.vector[0].pop(i)
        pointbox.checkWkbType()            
      return pointbox
    
    elif pointbox.category == 'Polygon':
      geom = None
      newVector = []
      for i in range(len(pointbox.vector)):
        lv2 = pointbox.vector[i]
        flag = False
        newRing = []
        for pt in lv2[0]:
          if math.fabs(pt.y()) < 0.0001:
            flag = True
            if pt.y() < 0:
              pt.setY(pt.y() + 0.001)
            else:
              pt.setY(pt.y() - 0.001)
          newRing.append(pt) 
      
        if flag:
          lv2[0] = newRing
          if geom == None:
            geom = dummy.fromPolygon(lv2)
          else:
            geom = geom.combine(dummy.fromPolygon(lv2))
        else:
          newVector.append(lv2)
      
      if geom != None:
        newPB = pointBox(geom)
        for lv2 in newVector:
          newPB.vector.append(lv2)
        newPB.checkWkbType()  
        return newPB
      else:
        return pointbox         
  
  '''
  Auxillary methods
  '''
     
  '''Counts the vertices in the given QgsGeometry (qgsGeom) '''
  def countVertices(self, qgsGeom):
    pb = pointBox(qgsGeom)
    count = 0
    for lv2 in pb.vector:
      for lv1 in lv2:
        count += len(lv1)
    return count

  '''
  Calculates the difference between the given longitude (lon) and
  the prime meridian. If the result overflowed from [-180, +180],
  it would be fixed to fit in the domain (eg. -270 --> +90)
  '''
  def normalisedLongitude(self, lon):
    lng = lon - self.lon0
    lng += 360 if lng < -180 else 0
    lng -= 360 if lng > 180 else 0
    return lng

  '''
  Transforms each point in the given pointBox to WGS84 lat/lon point
  '''
  def transformGeomToWGS84(self, pointbox):
    if self.epsg != 'EPSG:4326':
      for lv2 in pointbox.vector:
        for lv1 in lv2:
          for i in range(len(lv1)):
            pt = lv1[i]
            lv1[i] = self.xform.transform(pt)
    pointbox.checkWkbType()        

  '''
  Multiply the result (X, Y) of Gauss-Kruger equation by the scale 
  factor and then add the false easting to X and false northing to Y. 
  '''
  def applyScaleAndFalseNE(self, pointbox):
    for lv2 in pointbox.vector:
      for lv1 in lv2:
        for i in range(len(lv1)):
          pt = lv1[i]
          newX = pt.x() * self.scale + self.falseE
          newY = pt.y() * self.scale + self.falseN
          lv1[i] = QgsPoint(newX, newY)
  
  '''
  Reports "A point was transformed" to the progress bar sitting on the
  bottom of the user interface.
  '''
  def progressBar(self):
    if self.transformedPoints < self.vNum:
      self.transformedPoints += 1
      self.QProgBar.setValue(int(100 * self.transformedPoints / self.vNum))
    