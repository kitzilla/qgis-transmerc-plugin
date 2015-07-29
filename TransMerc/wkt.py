'''
*********************************************************************
wkt.py

Functions to extract projection parameters from OGC Well-Know Text(WKT)
encodings.

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

import re

'''A sample complex WKT for testing purposes'''
testWKT = 'PROJCS["WGS 84 / Antarctic Polar Stereographic",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137 ,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],PROJECTION["Polar_Stereographic"],PARAMETER["latitude_of_origin",-71],PARAMETER["central_meridian",0],PARAMETER["scale_factor",1],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","3031"]]'


'''
Extracts the parameters of map projections from strings of OGC's 
Well-Known Text (WKT) encodings.
'''
def crackWKT(wkt):
  rx1 = re.compile(r'([^\[\]]*)\[(.*)\]')
  rx2 = re.compile(r'"(.*?)",(.*)')
  rx3 = re.compile(r'authority', re.IGNORECASE)
  m = rx1.search(wkt)
  pwkt = None
  if m and len(m.groups())==2:
    n = rx2.search(m.group(2))
    if n and len(n.groups())==2:
      params = splitWKT(n.group(2))
      pwkt = wktProj(m.group(1), n.group(1))        
      for p in params:
        v = crackWKT(p)
        if v != None:
          pwkt.values.append(crackWKT(p))
    return pwkt 
  else:
    return wkt

'''
Splits a string which is a part of WKT by commas which are 
not enclosed by square brackets.

eg. INPUT : "[A,B],[C,D],[E,[F,G]]"
    OUTPUT: "[A,B]" | "[C,D]" | "[E,[F,G]]"
'''  
def splitWKT(string):
  res = [];  L = len(string)
  start = 0; braCnt = 0
  trim = re.compile(r'^\s+|\s+$')
  for i in range(L):
    c = string[i]
    braCnt += 1 if c == '[' else 0
    braCnt -= 1 if c == ']' else 0
    if c == ',' and braCnt == 0:
      if i > start:
        res.append(trim.sub('', string[start:i]))
      if i + 1 < len(string):
        start = i+1
    elif i == L-1:
      if L > start:
        res.append(trim.sub('', string[start:L]))
  return res    

'''
wktProj is a container box to store information extracted
by crackWKT()
'''
class wktProj():
  def __init__(self, ent, nam):
    self.entity = ent
    self.name = nam
    self.values = []
    self.indent = '&nbsp;&nbsp;'
    self.dindent = self.indent + '&nbsp;&nbsp;'

  def exportDescription(self):
    title = self.indent
    if self.entity.upper()=='PROJCS'  :
      title += 'Projected Coordinate System:<br>%s%s<br>' % (self.dindent, self.name)
    elif self.entity.upper()=='GEOGCS':
      title += 'Geographic Coordinate System:<br>%s%s' % (self.dindent, self.name)
    elif self.entity.upper()=='DATUM' :
      title += 'Datum: '+ self.name
    elif self.entity.upper()=='SPHEROID':
      title += 'Spheroid: %s<br>%sEquitorial Radius: %s [m]<br>%sFlattening: 1/%s' % (self.name, self.dindent, self.values[0], self.dindent, self.values[1])
    elif self.entity.upper()=='PRIMEM' :
      title += 'Prime Meridian: %s (%s)' % (self.values[0],self.name)
    elif self.entity.upper()=='UNIT' :
      title += 'Unit: '+ self.name
    elif self.entity.upper()=='PARAMETER' :
      title += self.name + ': '+ self.values[0]
    elif self.entity.upper()=='AXIS' :
      title += self.name + ': '+ self.values[0]
    else:
      return ''

    output = title + '<br>'
    if len(self.values) > 0:
      for v in self.values:
        if isinstance(v, wktProj):
          if len(v.indent) <= len(self.indent):
            v.indent = self.dindent
            v.dindent = v.indent + '&nbsp;&nbsp;' 
          output += v.exportDescription()
    return output      
