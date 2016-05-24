plotEQ<-function(Ldat, AQ, add=FALSE, prep=FALSE, TIT="UTM Projected Stations", proj=NULL, xlim=NULL, ylim=NULL )
  {
    zcols = c('red', 'blue', 'purple', 'cyan' )
    
    USTA = unique(Ldat$name)
    ww = !duplicated(Ldat$name)
    
    if(missing(proj)) {  proj = AQ$proj  }

    if(is.null(proj))  {  proj = AQ$proj  }
    
    EQ = AQ$EQ
####   get station X-Y values in km
    XY = GEOmap::GLOB.XY(Ldat$lat, Ldat$lon, proj)
    eXY = GEOmap::GLOB.XY(EQ$lat, EQ$lon, proj)

   
    if(!add)
      {
        plot(XY$x[ww],XY$y[ww], type='n' , asp=1, xlab="km", ylab="km", xlim=xlim, ylim=ylim )
      }

    if(prep)
      {
        return(NULL)
      }
    
    points( XY$x[ww],XY$y[ww] , bg="cyan", col='blue', pch=25 )
    points(eXY$x, eXY$y, col='green', pch=8)

    
    text(XY$x[ww],XY$y[ww], labels=Ldat$name[ww], pos=3, xpd=TRUE)

    qtip = ripper(AQ)

    trip = GEOmap::GLOB.XY(qtip[,1], qtip[,2], proj)
    
    points(trip$x, trip$y, col=zcols[qtip[, 5]], pch=qtip[, 5], cex=.5)
    lines(trip$x, trip$y,col=zcols[qtip[, 5]])

    title(main=TIT)
    
  }

