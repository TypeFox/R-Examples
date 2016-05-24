`GlobeView` <-
  function(phicen, lamcen, worldmap, MAXR, SEL=1, circol=rgb(1.0, .8, .8),  innercol="white",
           linecol=rgb(0, 0, 0),  mapcol=rgb(0, 0, 0), backcol="white" , add=FALSE )
  {
    
    if(missing(mapcol)) { eqcol=rgb(0, 0, 0) }
    if(missing(linecol)) { eqcol=rgb(0, 0, 0) }
    
    if(missing(circol)) { circol=rgb(1.0, .8, .8)  }
    if(missing(innercol)) { innercol=rgb(.8, .8, 1)  }

    
    if(missing(backcol)) { backcol="white" }
    if(missing(add)) { add=FALSE }

    
    if(missing(SEL)) { SEL = 1:length(worldmap$STROKES$index)  }

    phi1 =phicen
    lam0= lamcen
    mina = 20 
    maxa = MAXR+10
    DEGRAD = 0.017453293
    
    phi1 = DEGRAD * phi1
    lam0 = DEGRAD * lam0
    
    j = MAXR
    lat = DEGRAD * (phicen-j)
    lon = lam0
    xy = lamaz.eqarea(phi1,  lam0,lat, lon)
    r = sqrt(xy$x^2 + xy$y^2)
    ##
    ##C= circle()

    ##   i = pi * seq(from = 180, to = 360+180, by = 1)/180
    i = pi * seq(from = 0, to = 360, by = 1)/180

    
    C  = list(x=cos(i), y = sin(i))
    
    
    C$x = c(C$x, C$x[1])
    C$y = c(C$y, C$y[1])



    j=maxa
    
    
    lat = DEGRAD * (phicen-j)
    lon = lam0
    xy = lamaz.eqarea(phi1,  lam0,lat, lon)
    r = sqrt(xy$x^2 + xy$y^2)
    ## points(xy$x, xy$y, col=2)
    ##text(xy$x, xy$y, labels=j, pos=3)

    
    BDRx =  c(r*C$x)
    BDRy =  c(r*C$y)

    if(add==FALSE)
      {
        plot(range(BDRx)  , range( BDRy) , asp=TRUE, type='n', axes=FALSE, xlab='', ylab='', xaxs='i',yaxs='i',)
        polygon(BDRx, BDRy, col = innercol)
      }

##############  plot the circles
### for(j in seq(from=mina, to=maxa, by=20) )


###############   plot the map

    
    
    for(i in SEL)
      {
        
###  J = (worldmap$STROKES$index[i]+1):(worldmap$STROKES$index[i]+worldmap$STROKES$num[i])
        j1 = worldmap$STROKES$index[i]+1
        j2 = j1+worldmap$STROKES$num[i]-1
        
        lat = DEGRAD*worldmap$POINTS$lat[j1:j2]
        lon = DEGRAD*worldmap$POINTS$lon[j1:j2]
        
        
        xy = lamaz.eqarea(phi1,  lam0,lat, lon)
        
        are = sqrt(xy$x^2 + xy$y^2)
        
        flag = are>r
###  xy$x[flag] = NA
###  xy$y[flag] = NA
        
        if(!is.na(mapcol))     polygon(xy$x, xy$y, col=mapcol)
        
       if(!is.na(linecol))    lines(xy$x, xy$y, col=linecol)
###   numblats=worldmap$STROKES$num[i]
###    text( xy$x[numblats], xy$y[numblats] , labels=i, pos=3)
###    text( xy$x[1], xy$y[1], labels=i, pos=4)
###    cat(paste(sep=" ",i,",") )

        
        
      }



    antipolygon(BDRx, BDRy, col = backcol, corner=4)
    
    lines(BDRx, BDRy, col=circol , lwd=2, xpd=TRUE)
    
    ## grid.circle(0, 0, r)

    
    invisible(list(perx=BDRx, pery=BDRy))



    

    

############# 
    
  }

