parse.pde <-
function(card)
{
    ##### parse a pde card from the NEIC  website catalogue
    ####   format:
    #### PDE-W 2003   12 25 204233.72 -22.25  169.49  10 6.50 MwHRV
    #####c...|....1....|....2....|....3....|....4....|....5....|....6....|....7.$
  yr =  as.numeric(substr(card,9,12))
  mo = as.numeric(substr(card, 15,16))
  day  = as.numeric(substr(card,18, 19 ))
  hr  = as.numeric(substr(card, 21,22))
  mi  = as.numeric(substr(card, 23,24))
  sec  = as.numeric(substr(card, 25,29))
  lat= as.numeric(substr(card, 30, 36))
  lon= as.numeric(substr(card, 37, 44))
  depth= as.numeric(substr(card, 46, 48))
  
  if(is.numeric(depth))	{ z = depth } else { z = 0 }
  
  mag = as.numeric(substr(card, 50, 53))
  jd = getjul(yr, mo, day)
  locdate = list(yr=yr, jd=jd, mo=mo, dom=day, hr=hr, mi=mi, sec=sec,  lat=lat, lon=lon, depth=depth, z=z, mag=mag)
  return(locdate)
}

