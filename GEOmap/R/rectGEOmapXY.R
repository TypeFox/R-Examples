rectGEOmapXY<-function(lat=0, lon=0, PROJ=NULL, PMAT=NULL, ... )
  {
    if(missing(PMAT)) { PMAT=NULL }

    pxy = GLOB.XY(lat, lon, PROJ)


    if(!is.null(PMAT))
      {
        tem1 = trans3d(pxy$x, pxy$y, rep(0, length(pxy$y)) , PMAT)
      }
    else
      {
        tem1 =pxy

      }


    for(i in seq(from=1, to=length(tem1$x)-1, by=2))
      {
        rect(tem1$x[i],tem1$y[i],tem1$x[i+1],tem1$y[i+1], ...)
      }
    


  }

