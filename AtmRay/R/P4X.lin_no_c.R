P4X.lin_no_c = function(x, y, zs, zr, ATM = CheckAtm.lin(), maxerror = 3){
  # This function finds the ray parameter and azimuth that make the ray strike closest to the position (x, y) using Min2D.
  # x, y: horizontal position of target
  # c: intrinsic sound speed (vector, length n)
  # wx: east wind velocity (vector, length n)
  # wy: north wind velocity (vector, length n)
  # z: layer boundary heights (vector, length n+1)
  # zlim: c(ray_bottom, ray_top)
#  u0 = sin(az*pi/180) * ATM$wx0 + cos(az*pi/180) * ATM$wy0
#  v0 = sin(az*pi/180) * ATM$wy0 - cos(az*pi/180) * ATM$wx0
#  c0_eff = ATM$c0 + u0
#  gu = sin(az*pi/180) * ATM$gwx + cos(az*pi/180) * ATM$gwy
#  gv = sin(az*pi/180) * ATM$gwy - cos(az*pi/180) * ATM$gwx
#  gc_eff = ATM$gc + gu

  # determine p, azimuth range to search over

# This section determines minima of c_eff over the [zr, zs] range.  
# The next 6 lines define convenient, temporary variables for the calculation.
#  A = ATM$wx0*ATM$gwx + ATM$wy0*ATM$gwy
#  B = ATM$gwx^2 + ATM$gwy^2
#  C = ATM$wx0^2 + ATM$wy0^2
  
#  X0 = A^2 - ATM$gc^2 * C
#  X1 = 2*A*B - 2*ATM$gc^2*A
#  X2 = B^2 - ATM$gc^2 * B
   
#  zd = unique((-X1 + c(-1,1) * sqrt(X1^2 - 4*X0*X2))/(2*X2)) # points where dceff/dz = 0

#  zcrit = sort(c(zr, zs)) # elevation range in which the ray propagates
#  zcrit = c(zcrit, zd[zd > zcrit[1] & zd < zcrit[2]]) # add zd to this if it's in the range

#  cmin = min(ATM$c0 + ATM$gc * zcrit - sqrt((ATM$wx0 + ATM$gwx * zcrit)^2 + (ATM$wy0 + ATM$gwy * zcrit)^2))

#  pmax = 1/cmin # ray parameters greater than this cannot go from zs to zr, so don't test them
  l = 10
  thvec = seq(10^-3, 90 - 10^-3, length.out = l)
  azvec = seq(0, 360, length.out = l+1)[1:l]
  
  # loop until error is small
  error = Inf
  n = 0
  while(min(error, na.rm = TRUE) > maxerror){ # 1 m is small enough
    n = n+1
    if(n>20){ # it very rarely (if ever) seems to take more than 4 when it will find an answer, so 6 should be safe.
      k = NaN
      azk = NaN
      break
    }
    # define the meshgrids for this iteration
    thmat = meshgrid(thvec,azvec)$x
    azmat = meshgrid(thvec,azvec)$y
    c0_eff = ATM$c0 + sin(azmat*pi/180) * ATM$wx0 + cos(azmat*pi/180) * ATM$wy0
    gc_eff = ATM$gc + sin(azmat*pi/180) * ATM$gwx + cos(azmat*pi/180) * ATM$gwy
    pmat = sin(thmat * pi/180)/(c0_eff + gc_eff * zr)

    # do the propagation/error calculations
    L = Prop.lin_no_c(pmat, azmat, zs, zr, ATM)
    error = sqrt((L$x - x)^2 + (L$y - y)^2)

    # locate spot with minimum error and identify p and az for it
    k = which.min(error)
    thk = thmat[k]
    azk = azmat[k]

    # define new interval (2 steps in each direction from current position)
    dth = thvec[3]-thvec[2] # 3 and 2 to avoid epsilon effects at the edges
    daz = azvec[2]-azvec[1]
    thvec = seq(thk-1*dth, thk+1*dth, length.out = l)
    azvec = seq(azk-1*daz, azk+1*daz, length.out = l)
  }
#  print(n)
  return(list(p = pmat[k], az = azk, error = min(error, na.rm = TRUE)))
}
    
    
