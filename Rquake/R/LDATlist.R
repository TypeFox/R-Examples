LDATlist<-function(g1, w1)
{
sec = g1$STAS$sec[w1]

N = length(sec)
Ldat =    list(
  name = g1$STAS$name[w1],
  sec = g1$STAS$sec[w1],
  phase = g1$STAS$phase[w1],
  lat=g1$STAS$lat[w1],
  lon = g1$STAS$lon[w1],
  z = g1$STAS$z[w1],
  err= g1$STAS$err[w1],
  yr = rep(g1$LOC$yr , times=N),
  jd = rep(g1$LOC$jd, times=N),
  mo = rep(g1$LOC$mo, times=N),
  dom = rep(g1$LOC$dom, times=N),
  hr =rep( g1$LOC$hr, times=N),
  mi = rep(g1$LOC$mi, times=N) )

return(Ldat)


}

