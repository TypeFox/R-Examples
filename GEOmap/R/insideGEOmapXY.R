insideGEOmapXY<-function(lat, lon, PROJ=NULL, R=NULL, PMAT=NULL)
  {
    if(missing(PMAT)) { PMAT=NULL }
    if(missing(R)) { R = NULL }
####################   R =list(lat, lon, radius)
    
    pxy = GLOB.XY(lat, lon, PROJ)
    if(!is.null(PMAT))
      {
        tem1 = trans3d(pxy$x, pxy$y, rep(0, length(pxy$y)) , PMAT)
      }
    else
      {
        tem1 =pxy

      }

    u = par('usr')

    if(!is.null(R))
      {
        Rxy = GLOB.XY(R$lat, R$lon, PROJ)
        w = which( sqrt((Rxy$x-tem1$x)^2+(Rxy$y-tem1$y)^2)< R$radius )

        
      }
    else
      {
        w = which( tem1$x>u[1] & tem1$x<u[2] & tem1$y>u[3] & tem1$y<u[4])
      }

    
    return(w)
    
  }
