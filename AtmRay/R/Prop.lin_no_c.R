Prop.lin_no_c = function(p, az, zs, zr, ATM = CheckAtm.lin(list())){

  # Check Atmosphere and adjust heights to profile
  ATM = CheckAtm.lin(ATM)
  zs = zs - ATM$z0
  zr = zr - ATM$z0
  
  # Define effective sound speed slope and intercept
  u0 = cos(az*pi/180) * ATM$wy0 + sin(az*pi/180) * ATM$wx0
  v0 = sin(az*pi/180) * ATM$wy0 - cos(az*pi/180) * ATM$wx0
  c0_eff = ATM$c0 + u0
  gu = cos(az*pi/180) * ATM$gwy + sin(az*pi/180) * ATM$gwx
  gv = sin(az*pi/180) * ATM$gwy - cos(az*pi/180) * ATM$gwx
  gc_eff = ATM$gc + gu

  # Define source/receiver effective sound speed
  cs = c0_eff + zs * gc_eff
  cr = c0_eff + zr * gc_eff
#  c = c(cr, cs) # not necessary
  # These are useful variables to define separately
  epss = sqrt(1-p^2*cs^2)
  epsr = sqrt(1-p^2*cr^2)
  
  # Calculate displacement, time, amplitude
  # Equations from Jeff's notes for Villarrica paper from Slotnick
  rad = (zs - zr) * (epss - epsr)/(cr - cs)/p 
  T = abs((zs - zr)/(cs - cr) * log((cs * (1+epsr))/(cr * (1+epss))))

  dpdx = (cs - cr)*p^2*epsr*epss/ ((zs - zr) * (epss^2*epsr - epsr^2*epss + epsr*p^2*cs^2 - epss*p^2*cr^2))
  E = p*cs^2/(rad * (1-p^2*cs^2)) * dpdx
  # if density isn't defined in ATM, define it

  lambdar = (ATM$c0 + ATM$gc * zr)^2 * (ATM$rho0 + ATM$grho * zr) # the lame' parameter at the receiver
  A = sqrt(lambdar * E)
  tan = gv/gc_eff^2 * asin(cs*p)/p + (c0_eff*gv - v0*gc_eff)/gc_eff^2 * log((1/p + sqrt(1/p^2 - cs^2))/cs) - (gv/gc_eff^2 * asin(cr*p)/p + (c0_eff*gv - v0*gc_eff)/gc_eff^2 * log((1/p + sqrt(1/p^2 - cr^2))/cr))  # did this by hand over multiple days--int v dT = int v dT/dz dz
  x = sin(az*pi/180) * rad - cos(az*pi/180) * tan
  y = sin(az*pi/180) * tan + cos(az*pi/180) * rad
 # print(c(rad,tan))
  
  return(list(x = x, y = y, t = T, A = A, p = p))
}




#  here's what i did by hand.

#  z = 0:100
 # c0_eff = 333
#  gc_eff = 0.01
#  v0 = 5
#  gv = 0.1
#  c = c0_eff + gc_eff*z
 # p = 0.002
  

#  X2 = p*gv/gc_eff^2 * (-c*sqrt(1/p^2 - c^2)/2 + asin(c*p)/(2*p^2))

#  X1 = p*(c0_eff*gv-v0*gc_eff)/gc_eff^2 * sqrt(1/p^2 - c^2)

#  Z0 = gv*p/gc_eff^2 * (c*sqrt(1/p^2 - c^2)/2 + asin(c*p)/(2*p^2)) # note canceling terms in X2 and Z0

#  Z1 = p*(c0_eff*gv - v0*gc_eff)/gc_eff^2 * (-sqrt(1/p^2 - c^2) + log((1/p + sqrt(1/p^2 - c^2))/c)/p) # note canceling terms in Z1 and X1
  
#  X1 + Z1 + X2 + Z0
  
#  ATM = list(z=z,c=c[1:100],wy=rep(0, 100), wx = v0 + gv*z[1:100])
#  Prop.arb(0.002, 0, 100, 0, ATM)

#  z = seq(0, 100, 0.01)
#  n = length(z)
#  c0_eff = 333
#  gc_eff = 0.01
#  v0 = 5
#  gv = 0.1
#  c = c0_eff + gc_eff*z
#  p = 0.002
#  ATM = list(z=z,c=c[1:(n-1)],wy=rep(0, (n-1)), wx = v0 + gv*z[1:(n-1)])
#  Prop.arb(0.002, 0, 100, 0, ATM)$end
  
#  sum(diff(gv/gc_eff^2 * asin(c*p)/p + (c0_eff*gv - v0*gc_eff)/gc_eff^2 * log((1/p + sqrt(1/p^2 - c^2))/c)))

  
  
