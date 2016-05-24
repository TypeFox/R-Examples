`unpackAcard` <-
function(AC)
{
  ### unpackAcard(AC)
#### A 200508040604 45.19 64N0348  21W1390  1.91      12/007 116  0 0.02  0.1AB E1"
#### A 200406270357 30.03  0N4156  77W5069 12.00  4.1  0/000   0  0 0.00  0.0XX EC
#### c...|....1....|....2....|....3....|....4....|....5....|....6....|....7.$
  yr = as.numeric(substr(AC, 3, 6))
  mo = as.numeric(substr(AC, 7, 8))
  dom = as.numeric(substr(AC, 9, 10))
  hr = as.numeric(substr(AC, 11, 12))
  mi = as.numeric(substr(AC, 13, 14))
  sec  = as.numeric(substr(AC, 15, 20))
  LAT1 = as.numeric(substr(AC, 21, 23))
  LATNS = substr(AC, 24, 24)
  LAT2 = as.numeric(substr(AC, 25, 28))
  lat = LAT1+LAT2/(6000)


  jd = getjul(yr, mo, dom)


  LON1 = as.numeric(substr(AC, 30, 32))
  LONEW = substr(AC, 33, 33)
  LON2 = as.numeric(substr(AC, 34, 37))
  lon = LON1+LON2/(6000)

lat[LATNS=="S"] = -lat[LATNS=="S"]
lon[LONEW=="W"] = -lon[LONEW=="W"]
 
MAG  = substr(AC, 45, 48)

MAG = as.numeric(MAG)

  Z =   as.numeric(substr(AC, 38, 43))

  
      gap  = as.numeric(substr(AC, 57, 59))
  if(!is.numeric(gap)) gap = 0
      delta  = as.numeric(substr(AC, 60, 62))
    if(!is.numeric(delta)) delta = 0
  rms =as.numeric(substr(AC, 63, 67))
  if(!is.numeric(rms)) rms = 0
   hozerr  =as.numeric(substr(AC, 68, 72))
  if(!is.numeric(hozerr)) hozerr = 0
  
  return(list(yr=yr, mo=mo, dom=dom, hr=hr, mi=mi, sec=sec, jd=jd, lat=lat, lon=lon, z=Z, mag=MAG,gap=gap,delta=delta,rms=rms,hozerr=hozerr   ))
}

