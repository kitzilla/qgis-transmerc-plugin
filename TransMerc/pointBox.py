'''
*********************************************************************
pointBox.py

pointBox class for vertex processing.

*********************************************************************
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

   Begin     : 2012-07-01
   Copyright : (C) 2012 Go Sato 
   Email     : go.sato@zaq1.net

*********************************************************************
pointBox is a class to facilitate iterative processes on constituent 
points of QgsGeometry class.

All types of geometries in Quantum GIS API are represented as a
three dimensional list 'vector' in this class.

By using this class, points in QgsGeometry object can be processed 
regardless of what geometry type (Point, Polyline, Polygon, MultiPoint, 
MultiPolyline or MultiPolygon) the object is, in the following fashion
  
  for Level2 in pointBox.vector  #Level2: Polygon or MultiPolyline
    for Level1 in Level2         #Level1: Ring in Polygon, Polyline or MultiPoint
      for Point in Level1
        theProcessYouWantToDo(Point)  
  # len(vector) is always 1 unless the geometry is MultiPolygon type.
  # len(Level2) is always 1 when the geometry is Polyline or Point types
  # len(Level1) is    1 for Point type
  #                >= 2 for MultiPoint, Polyline and MultiPolyline types
  #                >= 4 for Polygon and MultiPolygon types
  # the first and last items of Level1 are identical for Polygon and MultiPolygon types. 
'''

from qgis.core import *


class pointBox:
  '''
  Constructor
  Constructs the three-dimensional list 'vector' (see the commets on
  the top) from the given QgsGeometry object (qgsGeom)
  '''
  def __init__(self, qgsGeom):
    self.vector = []
    self.wkbType = qgsGeom.wkbType()
    self.category = ''
    self.factory = QgsGeometry()
    if self.wkbType == QGis.WKBMultiPolygon:
      self.vector = qgsGeom.asMultiPolygon()
      self.category = 'Polygon'
    elif self.wkbType == QGis.WKBPolygon:
      level2 = qgsGeom.asPolygon()
      self.vector = [level2]
      self.category = 'Polygon'
    elif self.wkbType == QGis.WKBMultiLineString:
      level2 = qgsGeom.asMultiPolyline()
      self.vector = [level2]
      self.category = 'Polyline'
    elif self.wkbType == QGis.WKBLineString:
      level1  = qgsGeom.asPolyline()
      self.vector = [[level1]]
      self.category = 'Polyline'
    elif self.wkbType == QGis.WKBMultiPoint:
      level1 = qgsGeom.asMultiPoint()
      self.vector = [[level1]]
      self.category = 'Point'
    elif self.wkbType == QGis.WKBPoint:
      level0 = qgsGeom.asPoint()
      self.vector = [[[level0]]]
      self.category = 'Point'
  
  '''
  Check and update the geometry type (self.wkbType) of the pointBox.
  Geometry could easily change from multi- to non multi- or vice virsa
  by spatial operations (eg. union, difference, intersection etc.) 
  It is recommended to call method each time after operations that
  may affect the topology of the geometry stored in the pointBox.
  '''
  def checkWkbType(self):
    allLv1AreClosed = True
    allLv1HaveMoreThan2Pts = True
    allLv1HaveMoreThan4Pts = True
 
    for lv2 in self.vector:  #lv2: Polygon or MultiPolyline
      for lv1 in lv2:        #lv1: Ring in Polygon, Polyline or MultiPoint
        if len(lv1) <= 1:
          allLv1HaveMoreThan2Pts = False
          allLv1HaveMoreThan4Pts = False
          allLv1AreClosed = False
          break        
        elif len(lv1) < 4:
          allLv1HaveMoreThan4Pts = False
          allLv1AreClosed = False
          break
        else:    
          pt1 = lv1[0];  pt2 = lv1[-1]
          if pt1 != pt2:
            allLv1AreClosed = False
            break

    if self.category == 'Polygon':
      if allLv1AreClosed and allLv1HaveMoreThan4Pts:
        self.wkbType = QGis.WKBMultiPolygon if len(self.vector) > 1 else QGis.WKBPolygon
      else:
        self.wkbType = QGis.WKBUnknown
    
    if self.category == 'Polyline':
      if allLv1HaveMoreThan2Pts:
        self.wkbType = QGis.WKBMultiLineString if len(self.vector[0]) > 1 else QGis.WKBLineString
      else:
        self.wkbType = QGis.WKBUnknown
    
    if self.category == 'Point':
      self.wkbType = QGis.WKBMultiPoint if len(self.vector[0][0]) > 1 else QGis.WKBPoint

  '''
  Append the vector of another pointBox to the vector.
  Merger between adjacent or overlapping geometries is not considered
  in this method and that is why this is named 'simply'
  '''
  def simplyUnite(self, pointbox):
    if self.category == pointbox.category:
      if self.category == 'Polygon':
        for lv2 in pointbox.vector:
          self.vector.append(lv2)
      elif self.category == 'Polyline':
        for lv1 in pointbox.vector[0]:
          self.vector[0].append(lv1)
      elif self.category == 'Point':
        for lv0 in pointbox.vector[0][0]:
          self.vector[0][0].append(lv0)    
        
  '''
  Export the vector to QgsGeometry object.
  Note that if there is any discrepancy between the reported geometry
  type (self.wkbType) and the actual structure of the vector, 
  you might not retain what you expect.
  To make self.wkbType in line with the vector structure, use checkWkbType()
  before export. 
  '''
  def toQGSGeometry(self):
    geom = 0
    if self.wkbType == QGis.WKBMultiPolygon:
      geom = self.factory.fromMultiPolygon(self.vector)
    elif self.wkbType == QGis.WKBPolygon:
      geom = self.factory.fromPolygon(self.vector[0])
    elif self.wkbType == QGis.WKBMultiLineString:
      geom = self.factory.fromMultiPolyline(self.vector[0])
    elif self.wkbType == QGis.WKBLineString:
      geom = self.factory.fromPolyline(self.vector[0][0])
    elif self.wkbType == QGis.WKBMultiPoint:
      geom = self.factory.fromMultiPoint(self.vector[0][0])
    elif self.wkbType == QGis.WKBPoint:
      geom = self.factory.fromPoint(self.vector[0][0][0])

    return geom



       