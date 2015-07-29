'''
*********************************************************************
__init__.py

Initialises TransMerc plugin
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
def name():
  return "TransMerc"
def description():
  return "Accurately and universally transforms the projection of vector layers to Transverse Mercator using non-approximated Gauss Kruger formulae."
def version():
  return "Version 0.1"
def qgisMinimumVersion():
  return "1.7"
def classFactory(iface):
  from gkPluginInterface import gkPluginInterface
  return gkPluginInterface(iface)
def authorName():
  return "Go Sato"

def icon():
    """Icon"""
    return "transmerc_icon1.png"