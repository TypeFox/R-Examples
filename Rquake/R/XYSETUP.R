XYSETUP <-
function(STAS, init, vel)
  {
    proj = GEOmap::setPROJ(type=2, LAT0 =median(STAS$lat) , LON0 = median(STAS$lon) )
    XY = GEOmap::GLOB.XY(STAS$lat, STAS$lon, proj)
    gloc = GEOmap::GLOB.XY(init[1], init[2], proj)

    h1 = c(gloc$x, gloc$y, init[3], init[4])
    
    XY$r = sqrt((XY$x-h1[1])^2 + (XY$y-h1[2])^2)

    XY$c = rep(0, length(XY$r))
    XY$s = rep(0, length(XY$r))

    
    XY$c[XY$r>0] = (XY$x[XY$r>0]-h1[1])/XY$r[XY$r>0]
    XY$s[XY$r>0] = (XY$y[XY$r>0]-h1[2])/XY$r[XY$r>0]

    N = length(STAS$phase)

    elcor =  rep(0, times=N)

    DZ = STAS$z - mean(STAS$z)

    elcor[STAS$phase=="P"] = DZ[STAS$phase=="P"]/vel$vp[1]
    elcor[STAS$phase=="S"] = DZ[STAS$phase=="S"]/vel$vs[1]

    XY$cor = elcor
    XY$phase = STAS$phase
    XY$sec = STAS$sec

    return(XY)


  }
