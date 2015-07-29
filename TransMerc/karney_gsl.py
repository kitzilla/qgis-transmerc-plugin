'''
*********************************************************************
karney_gsl.py

Gauss Kruger Functions for trasformation from lat/lon to XY coordinates 
of Transverse Mercator Projection.

*********************************************************************
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

   Begin     : 2012-07-01
   Copyright : (C) 2008-2012 Charles Karney
               (C) 2012 Go Sato 
   Email     : go.sato@zaq1.net

*********************************************************************
This code is ported from TransverseMercatorExact.cpp in GeographicLib
(http://geographiclib.sourceforge.net/) to Python. The original code
was created by of Charles Karney and licenced to port, modify and re-
distribute under MIT/X11 License.
This code itself is licensed under GNU-GPL License.

The following modifications was added to the original to make it sure 
to work as a part of a Quantum GIS plugin.

1) Jacobi elliptic functions were switched to link their equivalents 
   in GNU Scientific Library (GSL) DLL.
2) TransverseMercatorExact::Forward() in the original code was divided
   into two functions (Forward and intermediateUV) to verify this code
   by comparing with another Gauss-Kruger code based on (Dozier 1980)
   
References:
  Improved Algorithm for Calculation of UTM and Geodetic Coordinates
  Dozier, J., 1980, NOAA Technical Report NESS81
*********************************************************************  
'''


import math
import numpy
import ctypes

libgsl = ctypes.CDLL('gsl.dll')

# Complete Elliptic Integral of the 1st Kind from GSL
def ellipK(m):
  func = libgsl['gsl_sf_ellint_Kcomp']
  func.restype = ctypes.c_double  
  func.argtypes = (ctypes.c_double, ctypes.c_uint)
  return func(math.sqrt(m), 0)

# Complete Elliptic Integral of the 2nd Kind from GSL
def ellipe(m):
  func = libgsl['gsl_sf_ellint_Ecomp']
  func.restype = ctypes.c_double  
  func.argtypes = (ctypes.c_double, ctypes.c_uint)
  return func(math.sqrt(m), 0)

# Incomplete Elliptic Integral of the 2nd Kind from GSL
def ellipeinc(phi, m):
  func = libgsl['gsl_sf_ellint_E']
  func.restype = ctypes.c_double  
  func.argtypes = (ctypes.c_double, ctypes.c_double, ctypes.c_uint)
  return func(phi, math.sqrt(m), 0)
  
# Jacobian Elliptic Functions from GSL
def ellipj(u, m):
  func = libgsl['gsl_sf_elljac_e']
  func.restype = ctypes.c_int  
  func.argtypes = (ctypes.c_double, ctypes.c_double, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double))
  sn = ctypes.c_double()
  cn = ctypes.c_double()
  dn = ctypes.c_double()
  func(u, m, ctypes.byref(sn), ctypes.byref(cn), ctypes.byref(dn))
  return [sn.value, cn.value, dn.value, math.atan2(sn.value, cn.value)]
  

class TransverseMercatorExact:
  def __init__(self, equatorRadius, polarRadius, originLongitude, originScale, extendp):
    self.a = equatorRadius
    self.b = polarRadius
    self.mu = (self.a**2 - self.b**2) / self.a**2
    self.mv = 1.0 - self.mu
    self.f = math.fabs(self.a - self.b) / min(self.a, self.b)
    self.lon0 = originLongitude    
    self.k0 = originScale
    self.e  = math.sqrt(self.mu)             # equiv to k in plugin
    self.ep2 = self.mu / self.mv             # m / (1 - m)
    self.extendp = extendp
    self.Eu = ellipe(self.mu)                # complete elliptic integral (2nd) of mu
    self.Ev = ellipe(self.mv)                # complete elliptic integral (2nd) of mv
    self.Ku = ellipK(self.mu)
    self.Kv = ellipK(self.mv)
    self.tolerance = numpy.finfo(float).eps
    self.tolerance2 = 0.1 * self.tolerance;
    self.taytol    = self.tolerance ** 0.6
    self.overflow  = numpy.finfo(float).max
  
  def taup(self, tau):
    tau1 = math.hypot(1.0, tau)
    sig = math.sinh( self.e * math.atanh(self.e * tau / tau1))
    return math.hypot(1.0, sig) * tau - sig*tau1;
    
    
  def taupinv(self, taup):
    tau = taup / self.mv
    stol = self.tolerance * max(1.0, math.fabs(taup))
    # min iterations = 1, max iterations = 2; mean = 1.94
    for i in range(100):
      tau1 = math.hypot(1.0, tau)
      sig = math.sinh( self.e * math.atanh(self.e * tau / tau1 ))
      taupa = math.hypot(1.0, sig) * tau - sig * tau1
      dtau = (taup - taupa) * (1.0 + self.mv * tau**2.0) 
      dtau /= self.mv * tau1 * math.hypot(1.0, taupa)
      tau += dtau
      if math.fabs(dtau) < stol:
        break
    return tau
  
  def zeta(self, snu, cnu, dnu, snv, cnv, dnv):
    '''
     Lee 54.17 but write
     atanh(snu * dnv) = asinh(snu * dnv / sqrt(cnu^2 + _mv * snu^2 * snv^2))
     atanh(_e * snu / dnv) =
             asinh(_e * snu / sqrt(_mu * cnu^2 + _mv * cnv^2))
    '''
    d1 = math.sqrt(cnu*cnu + self.mv * (snu*snv)**2.0)
    d2 = math.sqrt(self.mu * cnu*cnu + self.mv * cnv*cnv)
    t1 = snu * dnv
    if d1 != 0.0:
      t1 /= d1
    else: 
      t1 = self.overflow if snu >= 0.0 else -self.overflow
    
    t2 = self.e * snu
    if d2 != 0.0:
      t2 = math.sinh(self.e * math.asinh(t2 / d2))
    else:
      t2 = self.overflow if snu >= 0.0 else -self.overflow
    
    taup = t1 * math.hypot(1.0, t2) - t2 * math.hypot(1.0, t1);
    lam = 0.0
    if d1 != 0.0 and d2 != 0.0:
      lam = math.atan2(dnu*snv, cnu*cnv) - self.e * math.atan2(self.e*cnu*snv, dnu*cnv)
    return [taup, lam]
  
  def dwdzeta(self, snu, cnu, dnu, snv, cnv, dnv):
    '''
     Lee 54.21 but write (1 - dnu^2 * snv^2) = (cnv^2 + _mu * snu^2 * snv^2)
     (see A+S 16.21.4)
    '''
    d = self.mv * (cnv**2.0 + self.mu * (snu*snv)**2.0)**2.0
    du =  cnu*dnu*dnv * (cnv**2.0 - self.mu * (snu*snv)**2.0) / d;
    dv = -snu*snv*cnv * ((dnu*dnv)**2.0 + self.mu * cnu**2.0) / d; 
    return [du, dv]
  
  def zetainv0(self, psi, lam):
    retval = False
    psi0 = self.e * math.pi * 0.5
    lam0 = (1.0 - 2.0*self.e) * math.pi * 0.5
    dlam = lam - (1.0 - self.e) * math.pi*0.5
    u = 0.0;  v = 0.0
    
    if psi < -psi0 * 0.5 and lam > lam0 and psi < dlam:
      '''
       N.B. this branch is normally not taken because psi < 0 is converted
       psi > 0 by Forward.
       There's a log singularity at w = w0 = Eu.K() + i * Ev.K(),
       corresponding to the south pole, where we have, approximately
         psi = _e + i * pi/2 - _e * atanh(cos(i * (w - w0)/(1 + _mu/2)))
       Inverting this gives:
      '''
      psix = 1.0 - psi / self.e
      lamx = (math.pi * 0.5 - lam) / self.e
      u = math.sin(lamx) / math.hypot(math.cos(lamx), math.sinh(psix))
      u = math.asinh(u) * (1.0 + 0.5*self.mu)
      v = math.atan2(math.cos(lamx), math.sinh(psix)) * (1.0 + 0.5*self.mu);
      u = self.Ku - u;
      v = self.Kv - v;
    elif psi < psi0 and lam > lam0:
      '''
       At w = w0 = i * Ev.K(), we have
           zeta = zeta0 = i * (1 - _e) * pi/2
           zeta' = zeta'' = 0
       including the next term in the Taylor series gives:
       zeta = zeta0 - (_mv * _e) / 3 * (w - w0)^3
       When inverting this, we map arg(w - w0) = [-90, 0] to
       arg(zeta - zeta0) = [-90, 180]
      '''
      rad = math.hypot(psi, dlam)
      '''
       atan2(dlam-psi, psi+dlam) + 45d gives arg(zeta - zeta0) in range
       [-135, 225).  Subtracting 180 (since multiplier is negative) makes
       range [-315, 45).  Multiplying by 1/3 (for cube root) gives range
       [-105, 15).  In particular the range [-90, 180] in zeta space maps
       to [-90, 0] in w space as required.
      '''
      ang = (math.atan2(dlam-psi, dlam+psi) - 0.75 * math.pi) / 3.0
      
      # Error using this guess is about 0.21 * (rad/e)^(5/3)
      retval = rad < self.e * self.taytol
      rad = (3.0*rad / (self.mv * self.e)) ** (1.0 / 3.0)
      u = rad * math.cos(ang);
      v = rad * math.sin(ang) + self.Kv;
    else:
      '''
       Use spherical TM, Lee 12.6 -- writing atanh(sin(lam) / cosh(psi)) =
       asinh(sin(lam) / hypot(cos(lam), sinh(psi))).  This takes care of the
       log singularity at zeta = Eu.K() (corresponding to the north pole)
      '''
      v = math.asinh(math.sin(lam) / math.hypot(math.cos(lam), math.sinh(psi)))
      u = math.atan2(math.sinh(psi), math.cos(lam))
      
      # But scale to put 90,0 on the right place
      u *= 2.0 * self.Ku / math.pi;
      v *= 2.0 * self.Ku / math.pi;
    return [u, v, retval]
  
  def zetainv(self, taup, lam):
    psi = math.asinh(taup)
    
    scal = 1.0 / math.hypot(1.0, taup);

    (u, v, flag) = self.zetainv0(psi, lam)
    if flag:
      return [u, v]

    stol2 = self.tolerance2 / max(psi, 1.0)**2.0
    # min iterations = 2, max iterations = 6; mean = 4.0
    trip = 0
    
    for i in range(100):
      (snu, cnu, dnu, ph) = ellipj(u, self.mu)
      (snv, cnv, dnv, ph) = ellipj(v, self.mv)
      (tau1, lam1) =  self.zeta(snu, cnu, dnu, snv, cnv, dnv);
      (du1, dv1) = self.dwdzeta(snu, cnu, dnu, snv, cnv, dnv)
      
      tau1 -= taup;  lam1 -= lam;
      tau1 *= scal;
      delu = tau1 * du1 - lam1 * dv1
      delv = tau1 * dv1 + lam1 * du1
      u -= delu;  v -= delv;
      if trip > 0:
        break;
      delw2 = delu*delu + delv*delv;
      if delw2 < stol2 :
        trip += 1
    return [u, v]
  
  def sigma(self, snu, cnu, dnu, phu, v, snv, cnv, dnv, phv):
    '''
     Lee 55.4 writing
     dnu^2 + dnv^2 - 1 = _mu * cnu^2 + _mv * cnv^2
    '''
    d = self.mu * cnu*cnu + self.mv * cnv*cnv;
    
    xi = ellipeinc(phu, self.mu) - self.mu * snu * cnu * dnu / d
    eta = v - ellipeinc(phv, self.mv) + self.mv * snv * cnv * dnv / d
    return [xi, eta]
  
  def dwdsigma(self, snu, cnu, dnu, snv, cnv, dnv):
    '''
     Reciprocal of 55.9: dw/ds = dn(w)^2/_mv, expanding complex dn(w) using
     A+S 16.21.4
    '''
    d = self.mv * (cnv*cnv + self.mu * (snu*snv)**2.0)**2.0;
    dnr = dnu * cnv * dnv
    dni = -self.mu * snu * cnu * snv
    du = (dnr*dnr - dni*dni) / d;
    dv = 2.0 * dnr * dni / d;
    return [du, dv]
  
  def sigmainv0(self, xi, eta):
    retval = False;
    KEv = self.Kv - self.Ev
    u = 0.0; v = 0.0
    if eta > 1.25 * KEv  or (xi < -0.25 * self.Eu and xi < eta - KEv):
      '''
      sigma as a simple pole at w = w0 = Eu.K() + i * Ev.K() and sigma is
      approximated by
      sigma = (Eu.E() + i * Ev.KE()) + 1/(w - w0)
      '''
      x = xi - self.Eu 
      y = eta - KEv
      r2 = x*x + y*y
      u = self.Ku + x / r2
      v = self.Kv - y / r2
    elif (eta > 0.75 * KEv and xi < 0.25 * self.Eu) or eta > KEv :
      '''
      At w = w0 = i * Ev.K(), we have
          sigma = sigma0 = i * Ev.KE()
          sigma' = sigma'' = 0
      including the next term in the Taylor series gives:
      sigma = sigma0 - _mv / 3 * (w - w0)^3
      When inverting this, we map arg(w - w0) = [-pi/2, -pi/6] to
      arg(sigma - sigma0) = [-pi/2, pi/2]
      mapping arg = [-pi/2, -pi/6] to [-pi/2, pi/2]
      '''
      deta = eta - KEv
      rad = math.hypot(xi, deta)
      '''
      Map the range [-90, 180] in sigma space to [-90, 0] in w space.  See
      discussion in zetainv0 on the cut for ang.
      '''
      ang = (math.atan2(deta - xi, deta + xi) - 0.75 * math.pi) / 3.0
      # Error using this guess is about 0.068 * rad^(5/3)
      retval = rad < 2.0 * self.taytol 
      rad = math.pow(3.0 * rad / self.mv, 1.0/3.0)
      u = rad * math.cos(ang);
      v = rad * math.sin(ang) + self.Kv
    else:
      # Else use w = sigma * Eu.K/Eu.E (which is correct in the limit _e -> 0)
      u = xi  * self.Ku / self.Eu
      v = eta * self.Ku / self.Eu 
    return [u, v, retval]
  
  def sigmainv(self, xi, eta):
    (u, v, flag) = self.sigmainv0(xi, eta)
    if flag:
      return [u, v]
    
    # min iterations = 2, max iterations = 7; mean = 3.9
    trip = 0
    for i in range(1000):
      (snu, cnu, dnu, ph) = ellipj(u, self.mu)
      (snv, cnv, dnv, ph) = ellipj(v, self.mv)
      
      (xi1, eta1) = self.sigma(u, snu, cnu, dnu, v, snv, cnv, dnv)
      (du1, dv1 ) = self.dwdsigma(snu, cnu, dnu, snv, cnv, dnv)
      xi1 -= xi;
      eta1 -= eta;
      delu = xi1 * du1 - eta1 * dv1
      delv = xi1 * dv1 + eta1 * du1
      u -= delu;
      v -= delv;
      if trip > 0:
        break;
      delw2 = delu*delu + delv*delv
      if delw2 < self.tolerance2:
        trip += 1
        
    return [u, v]
  
  def Scale(self, tau, lam, snu, cnu, dnu, snv, cnv, dnv):
    sec2 = 1.0 + tau*tau
    '''
    Lee 55.12 -- negated for our sign convention.  gamma gives the bearing
    (clockwise from true north) of grid north
    '''
    gamma = math.atan2(self.mv *snu*snv*cnv, cnu*dnu*dnv)
    '''
    Lee 55.13 with nu given by Lee 9.1 -- in sqrt change the numerator
    from
       (1 - snu^2 * dnv^2) to (_mv * snv^2 + cnu^2 * dnv^2)
    to maintain accuracy near phi = 90 and change the denomintor from
       (dnu^2 + dnv^2 - 1) to (_mu * cnu^2 + _mv * cnv^2)
    to maintain accuracy near phi = 0, lam = 90 * (1 - e).  Similarly
    rewrite sqrt term in 9.1 as
       _mv + _mu * c^2 instead of 1 - _mu * sin(phi)^2
    '''
    k =  (self.mv *snv*snv + (cnu*dnv)**2.0) 
    k /= (self.mu *cnu*cnu + self.mv *cnv*cnv)
    k = math.sqrt(k) * math.sqrt(self.mv + self.mu/sec2) * math.sqrt(sec2)
    return [gamma, k]
  
  def intermediateUV(self, lat, lng):
    # Avoid losing a bit of accuracy in lon (assuming lon0 is an integer)
    lon = lng
    if lon - self.lon0 > 180:
      lon -= self.lon0 + 360
    elif lon - self.lon0 <= -180:
      lon -= self.lon0 - 360
    else:
      lon -= self.lon0;
    '''
    Now lon in (-180, 180]
    Explicitly enforce the parity
    '''
    latsign = 1 if lat >= 0.0 else -1
    latsign *= int(not self.extendp)
    lonsign = 1 if lon >= 0.0 else -1
    lonsign *= int(not self.extendp)
    
    lon *= lonsign;
    lat *= latsign;
    backside = (not self.extendp) and lon > 90;
    if backside:
      lon = 180 - lon
      if lat==0:
        latsign = -1;
    phi = math.radians(lat)
    lam = math.radians(lon)
    tau = math.tan(phi)
    tau =  self.overflow if phi >= 0 and tau < 0 else tau
    tau = -self.overflow if phi < 0 and tau >= 0 else tau

    # u,v = coordinates for the Thompson TM, Lee 54
    if lat == 90.0 :
      u = self.Ku
      v = 0.0
    elif lat == 0.0 and lon == 90.0 * (1.0 - self.e) :
      u = 0.0
      v = self.Kv
    else:
      (u, v) = self.zetainv(self.taup(tau), lam);
    return [u, v, latsign, lonsign, backside]
  
    
  def Forward(self, lat, lon):
    (u, v, latsign, lonsign, backside) = self.intermediateUV(lat, lon)
    (snu, cnu, dnu, phu) = ellipj(u, self.mu)
    (snv, cnv, dnv, phv) = ellipj(v, self.mv)
         
    (xi, eta) = self.sigma(snu, cnu, dnu, phu, v, snv, cnv, dnv, phv)
    
    if backside:
      xi = 2.0 * self.Eu - xi
    
    y = xi  * self.a * self.k0 * latsign
    x = eta * self.a * self.k0 * lonsign

    # Recompute (tau, lam) from (u, v) to improve accuracy of Scale
    (tau, lam) = self.zeta(snu, cnu, dnu, snv, cnv, dnv)
    tau = self.taupinv(tau);
    
    (gamma, k) = self.Scale(tau, lam, snu, cnu, dnu, snv, cnv, dnv)
    gamma = math.degrees(gamma)
    if backside:
      gamma = 180 - gamma
    gamma *= latsign * lonsign
    k *= self.k0
    return [x, y, gamma, k]
