plotnicetix<-function(nex, nwhy , proj, tlen=0.1, fonts=c("serif", "plain" ), PMAT=NULL, PLOT=TRUE )
  {
    if(missing(PLOT)) PLOT=TRUE
    if(missing(fonts))
      {
        tf <- "serif"
        fi <- "plain"
        
        fonts = c(tf, fi)
      }
    if(missing(tlen)) tlen = 0.1
    if(missing(PMAT)) { PMAT= NULL }

    
    rex = range(nex)
    rwhy = range(nwhy)
    
    gloc = XY.GLOB(rex, rwhy, proj)
    
    G = getnicetix(gloc$lat, gloc$lon)

    llim = XY.GLOB(rex, rwhy, proj)

    pusr = par("usr")
    pin = par("pin")
    
    ticlengthY = tlen*(pusr[4]-pusr[3])/(pin[2])
    ticlengthX = tlen*(pusr[2]-pusr[1])/(pin[1])
    



    for(i in 1:length(G$LON$DD))
      {
        
        ticlon = GLOB.XY(llim$lat[1], G$LON$DD[i], proj)

#####  this tick mark is slightly off - need to recalculate

        ##  points(ticlon$x[1], ticlon$y[1], pch=2)

        x1 =ticlon$x[1]
        y1 =ticlon$y[1]

        if(!is.null(PMAT))
          {
            tem1 = trans3d(x1, y1, 0 , PMAT)
            x1 = tem1$x
            y1 = tem1$y
          }
      
        
        x2 = x1
        y2 = y1-ticlengthY
        
        
        if(PLOT) segments(x1, y1, x2, y2)
        ## text(ticlon$x[1], ticlon$y[1]-ticlength, labels=G$LON$lab[i], pos=1, xpd=TRUE)

        asec = round(G$LON$sec[i])
        amin = round(G$LON$min[i])

        string = paste(sep="", format(G$LON$deg[i]), "\\de", format(G$LON$min[i]), "\\fm", format(asec), "\\sd")
        if(asec==0)
          { 
            string = paste(sep="", format(G$LON$deg[i]), "\\de", format(G$LON$min[i]), "\\fm")
          }
        if(asec==0 & amin==0)
          { 
            string = paste(sep="", format(G$LON$deg[i]), "\\de")
          }

        
      if(PLOT)  text(x2, y2,
             labels=string    ,
             pos=1, xpd=TRUE, vfont=fonts )
        ##  string = paste(format(G$LON$deg[i]), "\\de", format(G$LON$min[i]), "\\fm" , format(G$LON$sec[i]) 
        ##  text(ticlon$x[1], ticlon$y[1]-ticlength, string, vfont=c(tf, fi))

        PLONS = list(x=x2, y=y2, labels=string)

        
      }


    for(i in 1:length(G$LAT$DD))
      {
        
        ticlat = GLOB.XY( G$LAT$DD[i], llim$lon[1], proj)

        ###  this tic is not really good enough
        
        ##  points(ticlat$x[1], ticlat$y[1], pch=2)


        x1 =ticlat$x[1]
        y1 =ticlat$y[1]

        if(!is.null(PMAT))
          {
            tem1 = trans3d(x1, y1, 0 , PMAT)
            x1 = tem1$x
            y1 = tem1$y
          }
      
        x2 = x1-ticlengthX
        y2 = y1
        
        
       if(PLOT) segments(x1, y1, x2, y2)
        ## text(ticlat$x[1], ticlat$y[1]-ticlength, labels=G$LAT$lab[i], pos=1, xpd=TRUE)

        asec = round(G$LAT$sec[i])
        amin = round(G$LAT$min[i])

        string = paste(sep="", format(G$LAT$deg[i]), "\\de", format(G$LAT$min[i]), "\\fm", format(asec), "\\sd")
        if(asec==0)
          { 
            string = paste(sep="", format(G$LAT$deg[i]), "\\de", format(G$LAT$min[i]), "\\fm")
          }
        if(asec==0 & amin==0)
          { 
            string = paste(sep="", format(G$LAT$deg[i]), "\\de")
          }

        
      if(PLOT)  text(x2, y2,
             labels=string    ,
             pos=2, xpd=TRUE, vfont=fonts )

        PLATS = list(x=x2, y=y2, labels=string)


        
        ##  string = paste(format(G$LAT$deg[i]), "\\de", format(G$LAT$min[i]), "\\fm" , format(G$LAT$sec[i]) 
        ##  text(ticlat$x[1], ticlat$y[1]-ticlength, string, vfont=c(tf, fi)) 
      }

    invisible(list(lon=PLONS, lat=PLATS))


  }
