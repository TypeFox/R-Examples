contPFarrivals <- function(PF, stas, proj=NULL, cont=TRUE, POINTS=TRUE, image=FALSE ,
                            col=RSEIS::tomo.colors(50), gcol="black",   phase="P", add=TRUE)
  {
    ###  given a pickfile and a station file, contour arrivals
    
    if(missing(proj)) proj = NULL
   if(missing(cont)) cont=TRUE
   if(missing(POINTS)) POINTS=TRUE
   if(missing(image)) image=FALSE
    if(missing(col))  col=RSEIS::tomo.colors(50)


    NEX = 50

    if(length(PF$STAS$sec)<1)
      {
        print("ERROR: not enough station observations, no contouring")
        return(NULL)
      }

    
    pstas  = PF$STAS
    
    if(missing( phase))
      {
        ##  if phase is not specified try P, S, G
        if(any(pstas$phase=="P")) { phase = "P" }
        
        if(!any(pstas$phase=="P"))
          {
            if(any(pstas$phase=="S")) { phase = "S" }
          }
        
        if(!any(pstas$phase=="P") &  !any(pstas$phase=="S" ) )
          {
             phase = "G"
          }
        
        
      }
    
    w1 = which.min(pstas$sec)

    if(length(w1)<1) { return(NULL) }
    
    wp = which(pstas$phase==phase &   !is.na(pstas$lat) &  !is.na(pstas$lon)   )

    if(length(wp)<2) { return(NULL) }
    
    stan = pstas$name[wp]
    arr = pstas$sec[wp] - pstas$sec[w1]

   #  print("contPFarrivals: arr")
   #  print(arr)
    
    #  doAmap(stas)
    if(is.null(proj))  proj = GEOmap::setPROJ(type=2, LAT0 =median(stas$lat) , LON0 = median(stas$lon) )

    #############    convert the LATLON of stations to X-Y
    
    XY = GEOmap::GLOB.XY(stas$lat, stas$lon, proj)
    
    if(!add)
      {
        plot(XY, pch=6, xlab="km", ylab="km" , cex=.6 )

      }
    
   
   
    msta = match(stan, stas$name)

    ########   remove any non-station match
    msta = msta[!is.na(msta)]
    
    if(length(msta)<3)
      {
        print("ERROR: not enough matching stations, no contouring")
        return(NULL)
      }


    xy = GEOmap::GLOB.XY(  pstas$lat[wp] , pstas$lon[wp],  proj)

                         
    ex = seq(from=min(xy$x), to=max(xy$x), length=NEX)
    why = seq(from=min(xy$y), to=max(xy$y), length=NEX)


     DF = cbind(x=xy$x , y=xy$y ,  z=arr)
   ##   zed  = interp(x=xy$x , y=xy$y, z=arr, ex, why, duplicate="mean" )
     zed  = MBA::mba.surf(DF, NEX, NEX, n = 1, m = 1, h = 8, extend=FALSE)$xyz.est
    ##  ex = zed[[1]]
    ##  why = zed[[2]]

    if(image) image(zed, col=col , add=TRUE)
    if(POINTS)
      {
        text(XY, labels=stas$name, pos=3, cex=.6)
        points(XY$x[msta] , XY$y[msta], col='red', cex=1.1)
      }

    
    if(cont)
      {
        contour(zed, add=TRUE,  col='blue' )
    ##    contour(mzed, add=TRUE, col='black' )
        
      }
    return(proj)
  }
