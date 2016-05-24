`locworld` <-
function(shiftlon=0, col="brown", n=2)
  {
    if(missing(shiftlon)) { shiftlon=0 }
    if(missing(col)) { col= "brown" }
    if(missing(n)) {  n=2 }
    
    if(n<=0)
      {
        L = locator(type='p', pch=3, col=col)
      }
    else
      {
        L = locator(n, type='p', pch=3, col=col)
      }
    if(n==2){
    rect(L$x[1], L$y[1], L$x[2], L$y[2], border=col)}


    lon = RPMG::fmod(L$x+shiftlon, 360)
    lat = L$y
    
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


     A = list(lon=lon, lat=L$y  , LON=RPMG::fmod(lon, 360), LAT=L$y, utmbox=utmbox , x=L$x, y=L$y, UTM0=list(lam=utmlam0, phi=utmphi0), shiftlon=shiftlon)
   
    
    ## print(A)
    return(A)
  }

