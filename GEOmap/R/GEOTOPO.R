`GEOTOPO` <-
  function(TOPO, PLOC, PROJ, calcol=NULL, nx=500, ny=500, nb = 4, mb = 4, hb = 8, PLOT=TRUE)
  {


    if(missing(calcol)) {  CCOL = settopocol(); calcol = CCOL$calcol }
    if(missing(nx)) { nx=500  }
     if(missing(ny)) { ny=500  }
    
    if(missing(PLOT)) { PLOT=TRUE }

    if(is.null(calcol)) {  CCOL = settopocol(); calcol = CCOL$calcol  }



    cat("Extracting from Data Base....please wait....", file="", sep="\n")
    
    if(TRUE)
      {
###  here make adjustments so the topo part is square, not curved
###  to do this extend the borders of extraction

        nn = names(PLOC)
        
        ilon = grep("lon", nn, ignore.case = TRUE)
        ilat = grep("lat", nn, ignore.case = TRUE)
        
        if(length(ilon)<1)   { return(NULL) }
        if(length(ilat)<1)   { return(NULL) }

        
        A = list(lat=PLOC[[ilat[1]]], lon=PLOC[[ilon[1]]], LAT=PLOC[[ilat[1]]], LON=RPMG::fmod(PLOC[[ilon[1]]], 360) )
        
        
        
        PG  = GLOB.XY(A$lat, A$lon , PROJ  )

        ##   plot(PG, asp=1)

        
        dx = (PG$x[2]-PG$x[1])
        dy = (PG$y[2]-PG$y[1])

        pct = 10/100

        newLL = XY.GLOB(PG$x[1]-0.1*dx, PG$y[1]-pct*dy, PROJ  )
        newUR = XY.GLOB(PG$x[2]+0.1*dx, PG$y[2]+pct*dy, PROJ  )

        
        newPLOC = list(lat=c(newLL$lat,newUR$lat )  , lon=c(newLL$lon,newUR$lon )   )

        


        ZZ2 = subsetTOPO(TOPO, newPLOC, PROJ, nx=nx, ny=ny )
        
        
######   image(ZZ2)

        
        
        
      }
    
    


    if(PLOT==TRUE)
      {
        
        cat("Setting Colors....please wait....", file="", sep="\n")

       Cmat  = TOPOCOL(ZZ2$z, calcol)
        Dcol  = attr(Cmat, 'Dcol') 
        cat(".....plotting with persp....please wait....", file="", sep="\n")
        
        PMAT = persp(ZZ2$x, ZZ2$y, ZZ2$z, theta = 0, phi = 90, r=4000, col=Cmat[1:(Dcol[1]-1), 1:(Dcol[2]-1)] , scale = FALSE,
          ltheta = 120, lphi=30, shade = 0.75, border = NA, expand=0.001, box = FALSE )
        
      }
    else
      {
        PMAT = NA
      Cmat=NA
       Dcol=NA


      }
    
    
    invisible(list(PMAT=PMAT, xo=ZZ2$x, yo=ZZ2$y, IZ=ZZ2 ,Cmat=Cmat, Dcol=Dcol))
  }

