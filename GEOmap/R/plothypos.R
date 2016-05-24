plothypos<-function(lat, lon, z, proj=list(),  mag=NULL, cex=.4, pch =21, PMAT=NULL, alpha=NULL )
  {
    if(missing(PMAT)) { PMAT=NULL }
    if(missing(cex)) { cex=.4 }
    if(missing(alpha)) { alpha=NULL }
    if(missing(mag)) { mag=NULL }
     if(missing(pch)) { pch=21 }
    if(missing(proj)) {
      proj =  setPROJ(type = 2, LAT0=median(lat),LON0=median(lon) )

      }

    
    ####  use alpha = 0.3  

    kcol = c(0, 35, 70, 150, 300, 500, 800)
    
    r = c(255,  255,   66,   58,  255,  174)
    g = c(108,  163,  255,   75,   52,    0)
    b = c(86,   47,   82,  255,  253,    3)

    thecols=  rgb(r/255, g/255, b/255)

    

    if( !is.null(alpha)  )
      {
        thecols= adjustcolor(thecols, alpha.f=alpha)
      }

    if(!is.null(mag))
      {
        
        cex = getmagsize(mag, minsize=1,   slope=1, minmag=4,  maxmag=9, style=1)


      }
    

    ncols = length(thecols)
    
    h = which(z>kcol[length(kcol)])
    if(length(h)>0)
      {
        pointsGEOmapXY(lat[h], lon[h], PROJ=proj, col='black'   , bg=thecols[ncols]   , pch=21, cex=cex, PMAT=PMAT)
        
      }
    
    ################  events are plotted in order of color
    
    for(i in (length(kcol)):2 )
      {
        h = which(z>kcol[i-1] & z<=kcol[i])
        pointsGEOmapXY(lat[h], lon[h], PROJ=proj, col='black'   , bg=thecols[i-1]   , pch=21, cex=cex, PMAT=PMAT)
      }
    
  }




