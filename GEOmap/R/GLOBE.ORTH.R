GLOBE.ORTH<-function(lam0, phi1, R=1, plotmap=TRUE, plotline=TRUE, add=FALSE,
                     map='coastmap', mapcol =grey(.2) , linecol = grey(.7), fill=FALSE )
  {
    ###   use orthogonal projection to plot a global view
    if(missing(R)) R = 1
    if(missing(add)) add=FALSE
    if(missing( plotmap))  plotmap=TRUE
    if(missing( plotline))  plotline=TRUE

    
    if(missing(map))
      {
        plotmap = FALSE
   
      }

   if(missing( mapcol)) mapcol =grey(.2) 
   if(missing(linecol )) linecol = grey(.7) 
  
    
    Marids =  list()
    Parals =  list()


    lats = seq(from=90, to=-90, by=-2)
    lons = seq(from=0, to=340, by=20)

    for(i in 1:length(lons) )
      {
        Marids[[i]] = ortho.proj(lats, lons[i], lam0, phi1, R)
      }


    lats = seq(from=90, to=-90, by=-10)
    lons = seq(from=lam0-180, to=lam0+180, by=2)
    for(i in 1:length(lats) )
      {
        Parals[[i]] = ortho.proj(lats[i], lons , lam0, phi1, R)
      }

    ATOL= 0.4



    if(!add)
      {
        plot(c(-R, R), c(-R,R), asp=1, type='n', ann=FALSE, axes=FALSE)
      }

    
###   mapcol =rgb(1, .7, .7)
    if(plotmap)
      {

    for(i in 1:length(map$STROKES$num))
      {

        j1 = map$STROKES$index[i] + 1
        j2 = j1 + map$STROKES$num[i] - 1
        LONS = map$POINTS$lon[j1:j2]
        LATS = map$POINTS$lat[j1:j2]
        
        XY = ortho.proj(LATS, LONS , lam0, phi1, R)

        x  = XY$x[XY$cosc>=0]
        y = XY$y[XY$cosc>=0]

        n = length(x)

        if(n<=2) next
        
        dis = sqrt(  (x[1:(n-1)]-x[2:(n)])^2 + (y[1:(n-1)]-y[2:(n)])^2 )

        if(any(dis>ATOL))
          {

            w = which(dis>ATOL)
            x = insertNA(x, w )
            y = insertNA(y, w )
            ##  lines(nx, ny , col='green' )
          }
        lines(x, y , col=mapcol )
        if(fill==TRUE)
          {
            polygon(x,y,col=mapcol )

          }
      }
  }

    if(plotline)
      {

        circ = darc( rad=1, ang1=0, ang2=360, x1=0, y1=0, n=1)
        lines(circ)


        for(i in 1:length(Marids))
          {
            x = Marids[[i]]$x
            y = Marids[[i]]$y
            flag = Marids[[i]]$cosc>=0

            lines(x[flag],y[flag], col=linecol)


          }



        for(i in 1:length(Parals))
          {
            x = Parals[[i]]$x
            y = Parals[[i]]$y
            flag = Parals[[i]]$cosc>=0
            lines(x[flag],y[flag], col=linecol)
                                        # points(x[flag],y[flag], col='red')

                                        # locator()
            
          }

      }



  }
