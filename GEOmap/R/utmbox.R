`utmbox` <-
function(lat , lon)
  {
   
    utmx = seq(from=(-180), to=180, by=6)
    utmy = seq(from =-56, to =72, by=8)


    zlon = lon
     zlon[zlon>180 & zlon<=359.99] = zlon[zlon>180 & zlon<=359.99]-360

    
    htags =  LETTERS[6:23]
    htags = htags[-c(4, 10)]


    fzlon = findInterval(zlon[1], utmx)
    utmbox = list(x=fzlon, y = htags[findInterval(lat[1], utmy)], lon=utmx[fzlon] ,    lat=utmy[ findInterval(lat[1], utmy) ]   )

    fy = findInterval(lat[1], utmy)
    utmphi0 = (utmy[fy] +  utmy[fy+1])/2

    utmlam0 = (utmx[fzlon] +  utmx[fzlon+1])/2


     A = list(lon=lon, lat=lat  , LON=RPMG::fmod(lon, 360), LAT=lat, utmbox=utmbox , UTM0=list(lam=utmlam0, phi=utmphi0))
   
     ## utmbox(34.333333, 281)

    ## print(A)
    return(A)
  }

