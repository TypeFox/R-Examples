'getGEOperim'<-
  function(lon, lat, PROJ, N)
  {
    
    zim = GLOB.XY( rep(min(lat), length=N) ,seq(from=min(lon), to=max(lon), length=N),  PROJ )
    
    perim=zim
    
    
    zim = GLOB.XY( seq(from=min(lat), to=max(lat), length=N), rep(max(lon), length=N) ,  PROJ )
    perim=list(x=c(perim$x, zim$x), y=c(perim$y, zim$y))
    
    zim = GLOB.XY( rep(max(lat), length=N) ,seq(from=max(lon), to=min(lon), length=N),  PROJ )
    perim=list(x=c(perim$x, zim$x), y=c(perim$y, zim$y))
    
    zim = GLOB.XY( seq(from=max(lat), to=min(lat), length=N), rep(min(lon), length=N),  PROJ  )
    perim=list(x=c(perim$x, zim$x), y=c(perim$y, zim$y))
    
    perim=list(x=c(perim$x, perim$x[1]), y=c(perim$y, perim$y[1]))
    
    return(perim)
    
  }
