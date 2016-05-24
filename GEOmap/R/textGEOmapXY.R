textGEOmapXY<-function(lat=0, lon=0, labels=NULL, PROJ=NULL, PMAT=NULL, ... )
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


    text(tem1, labels=labels, ...)



  }

