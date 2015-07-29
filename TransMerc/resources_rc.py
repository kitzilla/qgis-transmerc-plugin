'''
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

   Begin     : 2012-07-01
   Copyright : (C) 2012 Go Sato 
   Email     : go.sato@zaq1.net
*********************************************************************
'''

from PyQt4 import QtCore

qt_resource_data = "\
\x00\x00\x06\x7e\
\x89\
\x50\x4e\x47\x0d\x0a\x1a\x0a\x00\x00\x00\x0d\x49\x48\x44\x52\x00\
\x00\x00\x16\x00\x00\x00\x16\x08\x02\x00\x00\x00\x4b\xd6\xfb\x6c\
\x00\x00\x00\x09\x70\x48\x59\x73\x00\x00\x0b\x13\x00\x00\x0b\x13\
\x01\x00\x9a\x9c\x18\x00\x00\x00\x04\x67\x41\x4d\x41\x00\x00\xb1\
\x8e\x7c\xfb\x51\x93\x00\x00\x00\x20\x63\x48\x52\x4d\x00\x00\x7a\
\x25\x00\x00\x80\x83\x00\x00\xf9\xff\x00\x00\x80\xe9\x00\x00\x75\
\x30\x00\x00\xea\x60\x00\x00\x3a\x98\x00\x00\x17\x6f\x92\x5f\xc5\
\x46\x00\x00\x05\xf4\x49\x44\x41\x54\x78\xda\x00\x43\x00\xbc\xff\
\x01\x0e\x6f\x35\x03\x1a\x10\xf3\xf0\xf1\xfc\x0e\xf1\x74\x18\x5e\
\x8c\x71\x7b\x13\x44\x34\x0e\x15\x0f\xff\xfd\xfe\xed\xe7\xee\x41\
\x56\x40\xb2\x94\xb0\x22\x32\x25\xe6\xe2\xe7\x3f\x2d\x33\xf5\x5f\
\x1a\xc4\x79\x88\x17\x35\x4c\xf7\xfa\xf2\x00\xff\xff\x01\x06\x03\
\xfd\xf0\xf9\x02\x00\x43\x00\xbc\xff\x01\x11\x89\x45\x0c\x6d\x3e\
\xe3\xc5\xb8\xd4\x44\xc4\x2b\x6a\xe3\xde\x31\xec\xbf\x9a\xb8\xaf\
\xcc\xab\x1c\x00\x1e\x0d\x00\x0b\x3a\x4f\x3f\xec\xe3\xeb\xdb\xce\
\xd8\xf6\x00\xf5\x21\x24\x24\x73\xdb\x8a\x58\x00\xba\xa9\xa1\x53\
\x1c\x3d\x6f\xf7\xf9\xf0\x01\x0e\x09\xfc\xcd\xe3\x02\x00\x43\x00\
\xbc\xff\x00\x03\x79\x35\x00\xbe\x44\x00\x79\x00\xd3\xff\xff\x00\
\x00\x00\x35\x00\x30\x7c\x50\x78\xd7\xcf\xd7\xaf\x97\xaf\x24\x00\
\x1b\x28\x00\x1f\xb6\xa3\xb6\x28\x00\x20\x05\x00\x00\xd8\xc0\xd8\
\x00\x00\x00\xff\xff\xff\x55\xff\xa9\x00\x85\x00\x19\xbe\x6b\x11\
\xc7\x64\x10\x9b\x4e\x02\x88\x81\xa1\xcb\x70\xe2\xff\x1b\xfb\xfe\
\xff\x00\xca\x85\x31\xb8\x45\xca\x04\x87\x0a\xf9\x59\x33\x18\x7a\
\xb3\xba\xfe\xf9\xf9\x57\x8a\x41\x00\x28\x5e\xec\x90\x07\x24\xeb\
\xfc\xaa\xfe\xbf\xfe\x0f\x33\x65\x2d\x43\x1b\x83\xd4\x85\x52\x81\
\xb9\x7e\x00\x01\xc4\x94\xc3\xe3\x9d\x57\x23\x76\x98\x2d\xbb\xd9\
\xba\x8c\x4d\x55\x98\x95\x99\x85\x5b\x80\x9b\x9d\x81\xeb\xe7\xef\
\x6f\x2c\xec\xcc\x4f\xff\xbd\x37\x65\x34\x78\x7c\xf9\xc1\xba\xbe\
\x55\xfd\x9b\x3a\x15\x44\x45\x81\xfa\x19\x19\x19\x57\x31\x04\x85\
\x56\xae\x7d\x76\x78\x22\x07\x1b\x3f\x40\x00\x31\x30\x08\x31\xdc\
\x98\x7c\x5c\x88\x81\xdb\x8b\xdd\x25\x5e\x21\x22\x52\xdc\xdf\x8d\
\xc1\xfe\xce\x85\x5b\xfb\x96\xef\x59\xdb\xb3\xfa\xe9\xcd\xc7\x65\
\xce\x05\x9f\x5e\x7f\x02\x5a\xbe\x6f\xe9\xde\x48\xb1\x00\x20\xe3\
\xea\x9e\xcb\xfa\x0c\x1a\x60\xb7\x6c\x61\xc8\x66\x00\x08\x20\x16\
\x61\x27\x93\xa3\x8b\xf6\xe7\x04\xe7\xdf\x38\x78\x95\x91\x89\xf1\
\xed\xcb\x0f\x7d\x47\xa7\xb8\x1b\xd8\x29\x33\xa8\x3f\x66\xb8\x5f\
\xf6\xb1\xee\xc0\xde\xbd\xfc\x79\xfc\x40\x0d\x1e\x4c\xf6\x67\xff\
\x9f\x39\xb9\xf9\xf8\x92\xae\xc5\x8b\x36\xaf\x00\xba\x05\x28\xe8\
\xdc\xd7\x03\x10\x00\x43\x00\xbc\xff\x01\x21\x69\x43\x2a\x97\xee\
\x8c\xcf\xa6\x64\x31\x5a\x07\x00\x0a\x9a\xd0\x9f\x48\x30\x40\x45\
\x2e\x48\x01\x03\x02\xc4\xcf\xc0\x5e\x60\x63\x02\x02\x03\x0a\x0e\
\x09\xf9\xf5\xfb\x6e\x9a\x71\x49\x04\x41\x04\x03\x06\xf9\xfa\xf5\
\x0c\x0e\x10\x76\x91\x7d\x39\x61\x38\xef\xe6\xeb\x02\x88\x05\x18\
\xff\xc0\xf8\xbb\x71\xe2\x86\x00\xab\x08\xd0\xff\xcf\xee\x3c\x63\
\x61\x60\x61\xe5\x60\x3b\x7c\x6a\x5f\x76\x49\x49\x5a\x77\xba\x38\
\xa3\x80\x09\x93\xd1\xb9\x7f\x67\x0f\xad\xdb\xb7\xfd\xcf\x7e\x56\
\x46\x46\x67\x06\xbb\x93\x9f\x8e\xb2\xc8\x30\x25\x32\x04\x28\x2a\
\x1a\x02\x04\x10\x03\xaf\xaf\x7e\x09\x43\x14\x30\xfe\x80\xe1\x0f\
\xf4\x67\x6f\x5c\xc7\xc1\x55\x07\xf8\x18\x18\xe7\x94\xcf\x72\x65\
\xb0\xd1\x62\x50\x39\xb2\xee\xb0\x0d\x83\xb1\x3d\x83\xe9\xfa\x09\
\x6b\x03\xf9\x3c\x80\x8e\x4f\xd7\x4d\x3c\xb3\xed\x74\x6d\x60\xb5\
\x29\x83\x92\x9e\xae\x35\x40\x00\x00\x43\x00\xbc\xff\x03\x48\x7d\
\x58\x4d\xfe\x35\xbf\xd9\xc5\x87\x84\x81\xf3\xff\xf8\x39\x51\x3c\
\x23\x2c\x21\xa5\x91\xab\x5b\x7b\x5d\xbd\xa3\xbc\xe2\xef\xde\x17\
\x0a\x18\x24\x2a\x25\xec\xe1\xee\xfe\xf9\xfc\xf7\xf5\xf9\xd5\xc5\
\xd3\x5a\x79\x5c\x33\x49\x35\xc1\xd7\xc2\x25\x07\x26\x0b\x42\x17\
\x02\x88\x81\xc1\x5c\x7e\x96\x51\xcb\xb6\x45\x5b\x41\xb1\xc5\xc0\
\x50\xe7\x5b\x19\x26\xea\xd7\x19\xd3\xaa\xcc\x20\xd5\x1a\xde\x00\
\x74\x91\x3c\x83\xb0\x3b\x83\xdd\x7f\x70\xc2\x93\x67\x10\xed\x8a\
\x6b\x73\x66\xb0\x05\x72\xc3\xa5\x43\x25\x19\x98\xe5\x44\xb4\x00\
\x02\x00\x43\x00\xbc\xff\x03\x21\x4d\x33\x18\xb3\xf9\x94\xaf\x99\
\xd5\xb1\xd3\x54\x44\x53\x66\x8a\x6b\xdc\xcf\xd8\x17\x20\x18\xeb\
\xc7\xe9\xe2\xd6\xe0\xf6\xf4\xf5\xfb\xdd\xfb\x75\xa3\x7d\x29\x2f\
\x23\xe3\xdc\xe4\x28\x1d\x2b\x48\x61\x4a\x6d\x95\x71\x01\x01\x01\
\x24\x00\x20\xcb\xda\xd3\xc7\xe6\xd1\x02\x88\x41\x2a\xc5\x36\x8d\
\xc1\x87\x15\x9c\xe6\x16\x35\x2c\x80\x24\xbe\x8e\x88\x66\x49\x06\
\x7e\xa0\xff\x03\xf8\x3d\x81\x5c\x36\x06\x86\x54\x8d\x18\x67\x06\
\x1b\x19\x06\xa1\xed\xb3\x36\x03\xd5\x3c\xbf\xfb\x0c\x18\x40\x2e\
\x0c\x26\x81\x31\x71\x00\x01\x00\x43\x00\xbc\xff\x04\xe7\xb2\xd8\
\x59\x77\x9b\x69\x44\x83\x7f\x4f\x7b\x90\xb2\x89\xbf\xa7\xbf\x3d\
\x01\x36\xc3\xff\xca\xcc\xbc\xcb\x4f\x6f\x53\xdc\xc7\xd8\x01\x01\
\x01\xe6\xdc\xe4\x5d\x24\x58\x1a\x1c\x17\x08\xe9\x09\x00\x03\x03\
\x2f\x00\x25\xd7\xc9\xd5\x06\x00\x08\x30\x1d\x2a\xc2\xda\xde\x02\
\x88\x21\xe7\xff\xfa\xff\x1e\x8f\x81\x06\xff\x7e\xfd\x57\x9e\x41\
\x04\x98\xfe\x44\x18\xb8\x80\xf1\x0f\xb4\xbc\x2d\xac\x11\x98\x9c\
\x14\x19\xc4\x57\x75\xaf\x00\xda\x0f\xca\x35\x7f\x7f\x2b\x33\x48\
\x78\x32\x39\xed\x59\xb2\xcb\x91\xc1\xaa\x50\x22\xd1\xde\xc9\x07\
\x20\x80\x58\xa6\x78\x07\x6a\x6f\xff\xff\x3f\xe6\x27\xa3\x28\x33\
\x24\xfd\x03\xd3\x5f\x73\x70\xbd\x00\x23\xdb\x87\xff\xbf\x6c\x57\
\x39\xee\x59\xb1\x43\xdf\x56\x77\xce\xac\x45\xff\x99\x99\xd4\x99\
\x65\x7d\x2d\x82\xf6\x9f\xd8\x2d\xa9\x24\xc9\xc4\xc0\xcc\xc2\xce\
\xca\xc2\xca\x0c\x10\x40\x0c\xe2\xad\x3e\x0c\xf3\x98\xa6\xfc\xff\
\xf9\x3f\xf8\x3d\xd0\x2d\xa7\xb7\x9e\xac\x72\x2f\x07\xe6\x54\x37\
\x06\xdb\x5a\xef\xf2\xd9\x85\xd3\x20\xe1\x0f\x94\x32\x60\x50\x05\
\xda\x0f\xf4\xc5\x9c\xb2\xe9\x91\x12\xc1\xa2\x0c\xbc\xc6\x0c\xaa\
\x1c\x4c\x62\x00\x01\xc4\x02\x2c\x7f\x64\xec\x2b\x73\x1a\xd9\x59\
\xd6\xfc\xfb\x2f\xb1\xa7\x6c\xfb\x26\x8f\x70\x17\x1e\x41\xae\xae\
\x15\x2d\xad\x59\x3d\xa6\xde\xe6\x10\xa7\x7d\x7c\xfc\xf9\xcc\xa1\
\x13\x40\x36\xd0\xfe\x96\xa0\x06\x0d\x5f\xdd\x65\xc9\x6b\x18\x1e\
\x32\x30\x68\x7d\x06\x08\x20\x06\xbe\x3a\x5b\xe1\x6b\xd9\xc0\xf2\
\x83\xa1\x91\x61\x16\x50\x89\xec\x51\xa0\x85\xdd\x11\x6d\xd7\xcf\
\x5c\xab\x09\xa8\x02\x15\x22\x92\x21\x11\xc2\xc1\xd6\x0c\x46\x4e\
\x0c\xd6\x40\xff\x03\x1d\x05\xb4\xff\xff\x89\xff\xd1\xff\x37\x9b\
\xfd\x9f\xce\x30\x99\x01\x20\x80\x58\x80\xe5\x1f\xb0\xfc\xfa\xf1\
\xf1\x83\x94\x7f\x69\xda\x04\xc6\xbf\x8f\xfe\xff\x4f\xfb\xcd\x3f\
\x9b\xe3\xfc\x8a\x43\xb7\x19\x6e\x39\xb7\xdb\x1c\x7e\xbe\xfe\xd9\
\x83\xaf\xae\x0a\x36\xc0\xf0\xff\xfb\xfb\x0f\xd0\xff\x1f\x1e\xbe\
\x00\xda\x7f\xdb\xfc\xc9\xa9\xed\xd3\x25\xc4\x12\x00\x02\x00\x43\
\x00\xbc\xff\x03\x05\x32\x19\x00\x05\x01\x00\xdf\xf1\xfb\xe8\xed\
\x10\x0d\x36\xa9\x93\x52\xf9\x9b\x4e\x48\xdc\x19\xb7\x00\xd4\xfd\
\x15\x0c\x0c\x43\x18\xd1\xc7\xe1\xc7\xe6\xd1\x05\x2c\x13\xd4\xfa\
\xe2\xec\x1b\xfa\xe9\x17\xf6\xcd\xdf\xdc\x7c\x7e\x68\x66\xaa\x3c\
\x18\x05\x49\xf7\xe0\xe2\x02\x0c\x00\x70\x33\x7e\x9c\x9c\x74\xcf\
\x80\x00\x00\x00\x00\x49\x45\x4e\x44\xae\x42\x60\x82\
"

qt_resource_name = "\
\x00\x13\
\x05\x26\x83\x27\
\x00\x74\
\x00\x72\x00\x61\x00\x6e\x00\x73\x00\x6d\x00\x65\x00\x72\x00\x63\x00\x5f\x00\x69\x00\x63\x00\x6f\x00\x6e\x00\x32\x00\x2e\x00\x70\
\x00\x6e\x00\x67\
"

qt_resource_struct = "\
\x00\x00\x00\x00\x00\x02\x00\x00\x00\x01\x00\x00\x00\x01\
\x00\x00\x00\x00\x00\x00\x00\x00\x00\x01\x00\x00\x00\x00\
"

def qInitResources():
    QtCore.qRegisterResourceData(0x01, qt_resource_struct, qt_resource_name, qt_resource_data)

def qCleanupResources():
    QtCore.qUnregisterResourceData(0x01, qt_resource_struct, qt_resource_name, qt_resource_data)

qInitResources()