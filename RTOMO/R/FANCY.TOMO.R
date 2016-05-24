`FANCY.TOMO` <-
function(MOD, i, COL=NULL, LIM=NULL, MAP=NULL, MAPLIM=NULL,
                     STA=NULL, staparams=list(col='green', pch=6, cex=.8, name=FALSE),
                     PTS=NULL, ptsparams=list(col='green', pch=6, cex=.8, name=FALSE),
                     TIT="Layer", mainTIT="Layer",  UNITS="", bkgr="DarkSlateGray4" )
  {

    

    if(missing(COL)) { COL=tomo.colors(100) }
    if(missing(LIM)) { LIM=NULL }
    if(missing(MAP)) { MAP=NULL }
    if(missing(MAPLIM)) { MAPLIM=NULL }
    
    if(missing(STA)) { STA=NULL }
    if(missing(staparams)) { staparams=list(col='green', pch=6, cex=.8, name=FALSE) }

    if(missing(PTS)) { PTS=NULL }
    if(missing(ptsparams)) { ptsparams=list(col='green', pch=6, cex=.8, name=FALSE) }

   
    if(missing(UNITS)) { UNITS="%" }
    if(missing(bkgr)) { bkgr="DarkSlateGray4" }
    
    if(missing(TIT)) { TIT= paste(sep=' ', "LAYER=", i, "Depth", MOD$D[i], "-", MOD$D[i+1]) }
    if(missing(mainTIT)) { mainTIT=NULL }

    if(is.null(LIM)) { LIM= range(MOD$MOD[[i]], na.rm=TRUE) }

   typefaces = c("serif", "sans serif", "script", "gothic english", 
        "serif symbol", "sans serif symbol")
    fontindeces = c("plain", "italic", "bold", "bold italic", 
        "cyrillic")
    typeface = typefaces[1]
    fontindex = fontindeces[1]
    vfont = c(typeface, fontindex)


    leestomo = TRUE
    if(leestomo)
      {
    
    DX = (MOD$x[2]-MOD$x[1])/2
    NX = length(MOD$x)
    
    PX = seq(from=MOD$x[1], to=MOD$x[NX]+((MOD$x[2]-MOD$x[1])) , length=NX+1)
    
    DY = (MOD$y[2]-MOD$y[1])/2
    NY = length(MOD$y)
    
    PY = seq(from=MOD$y[1], to=MOD$y[NY]+((MOD$y[2]-MOD$y[1])), length=NY+1)
    
  }
    
    pltomo(PX,PY,MOD$MOD,i, COL, zlim=LIM, bkgr=bkgr, xlab="km", ylab="km")

## GEOmap::ColorScale(MOD$MOD, loc=list(x=range(PX), y=range(PY)), offset=.2,
##          col = COL, units = "%Vel", font = 1, fontindex = 3, SIDE = 3)


##   GEOmap::ColorScale(MOD$MOD, loc=list(x=range(PX), y=range(PY)), offset=.5,
##            col = COL, units = "%Vel", font = 1, fontindex = 3, SIDE = 1)


    
   ##   HOZscale(LIM , col=COL  , units=UNITS, SIDE=1)
     
   
    
### HOZscale(MOD$MOD[[i]] , col=tomocolors  , units="%", SIDE=2)
        
        
###  image.SCALE(MOD$MOD[[i]] , col=rainbow(100), nlab=2)	
        ## PROJmap(JAPmap,  ADD=TRUE, COL=TRUE)
        if(!is.null(MAP))
          {
            ##  plotGEOmapXY(MAP, LIM=MAPLIM   , PROJ =MAP$PROJ ,   add=TRUE)
            GEOmap::plotGEOmapXY(MAP, PROJ=MAP$PROJ, add=TRUE, xpd=FALSE )
          }
        if(!is.null(STA))
          {
            if(is.null(STA$FLAG)) { sflag = rep(TRUE, length(STA$x))} else { sflag = STA$FLAG }
            
            GEOmap::pointsGEOmapXY(STA$lat[sflag], STA$lon[sflag], PROJ=MAP$PROJ, pch=staparams$pch, col=staparams$col, cex=staparams$cex)
            
            if(identical(staparams$name, TRUE))
              {
                
                GEOmap::textGEOmapXY(STA$lat[sflag], STA$lon[sflag], labels=STA$name[sflag], PROJ=MAP$PROJ, pos=3, cex=staparams$cex)
              }
            
            
            
          }
        
        if(!is.null(PTS))
          {
            if(is.null(PTS$FLAG)) { sflag = rep(TRUE, length(PTS$lat))} else { sflag = PTS$FLAG }
            GEOmap::pointsGEOmapXY(PTS$lat[sflag], PTS$lon[sflag], PROJ=MAP$PROJ, pch=ptsparams$pch, col=ptsparams$col, cex=ptsparams$cex)
            
            if(identical(ptsparams$name, TRUE))
              {
                
                GEOmap::textGEOmapXY(PTS$lat[sflag], PTS$lon[sflag], PROJ=MAP$PROJ, PTS$name[sflag], pos=3, cex=ptsparams$cex)
              }
        
          }
        
      
    
    
    ##	text(STA$x[STAFLAG], STA$y[STAFLAG], labels=STA$nam[STAFLAG], pos=3, col='green', cex=.65)
    ## print(TIT)


    loc=list(x=range(PX), y=range(PY))

    GEOmap::antipolygon(x=c(loc$x[1], loc$x[2], loc$x[2] , loc$x[1]),
                y=c(loc$y[1], loc$y[1], loc$y[2], loc$y[2]),col = "white" , corner = 1, pct = 0.4)

    rect(loc$x[1], loc$y[1], loc$x[2], loc$y[2])

    xax = pretty(range(PX))
    xax = xax[xax>=loc$x[1] & xax<=loc$x[2] ]

  yax = pretty(range(PY))
    yax = yax[yax>=loc$y[1] & yax<=loc$y[2] ]



    RSEIS::addtix(side=1, pos=loc$y[1], at=xax, tck = 0.8)
    text(xax, rep(loc$y[1],  length=length(xax)), labels=xax, pos=1, vfont=vfont)

    RSEIS::addtix(side=2, pos=loc$x[1], at=yax, tck = 0.8)
    text( rep(loc$x[1],  length=length(yax)), yax , labels=yax, pos=2, vfont=vfont)
   
    
  
     GEOmap::ColorScale(MOD$MOD, loc=list(x=range(PX), y=range(PY)), offset=.5,
          col = COL, units = "%Vel", font = 1, fontindex = 3, SIDE = 1)


    if(!is.null(TIT))
      {
        title(main=mainTIT, sub=TIT)
      }
    
    
  }

