Vlocate<-function(Ldat,EQ,vel, 
                  distwt = 10,
                  lambdareg =100,
                  REG = TRUE,
                  WTS = TRUE,
                  STOPPING = TRUE,
                  tolx = 0.1,
                  toly = 0.1,
                  tolz = 0.5,
                  RESMAX = c(.4,.5),
                  maxITER = c(7, 5, 7, 4),
                  PLOT=FALSE)
  {

###  Earthquake location program for local earthquakes.
    ##  this code includes a projection to cartesian coordinates.

    BADSOLUTION = FALSE
    Ksolutions = vector(mode="list")

#######   set up input list #########
####  set up the projection to local UTM coordinates
    MLAT = median(Ldat$lat)
    MLON = median(Ldat$lon)
    
    proj = GEOmap::setPROJ(type=2, LAT0=MLAT, LON0=MLON)

####   get station X-Y values in km
    XY = GEOmap::GLOB.XY(Ldat$lat, Ldat$lon, proj)
###   add to Ldat list
    Ldat$x = XY$x
    Ldat$y = XY$y
    

    
###########   Earthquake location starting guess
    if(is.na(EQ$lat))
      {
        cat("Guess at first arrival station", sep="\n" )
        
        wstart = which.min(Ldat$sec)
        EQ = list(lat=Ldat$lat[wstart], lon=Ldat$lon[wstart], z=6, t=Ldat$sec[wstart]-1)
        EQ$x = XY$x[wstart]
        EQ$y = XY$y[wstart] 
      }else{
        eqxy = GEOmap::GLOB.XY(EQ$lat, EQ$lon, proj)
        EQ$x = eqxy$x
        EQ$y = eqxy$y
      }


################### crude estimate  iterations  

###     /* fix depth for initial entry to locath */
  ###    maxITER = 7
###print(EQ)
    AQ = XYlocate(Ldat,EQ,vel, 
      maxITER = maxITER[1],
      distwt = distwt,
      lambdareg =lambdareg ,
      FIXZ =TRUE,
      REG = REG,
      WTS = WTS,
      STOPPING = STOPPING,
      RESMAX = c(0,0),
      tolx =   tolx,
      toly = toly ,
      tolz = tolz, PLOT=PLOT)

    Ksolutions[[1]] = AQ$guesses

    if(is.null(AQ$EQ))
       {
         return(list(EQ=NULL ))
 

       }
       
    eqLL = GEOmap::XY.GLOB(AQ$EQ$x, AQ$EQ$y, proj)
    AQ$EQ$lat = eqLL$lat
    AQ$EQ$lon = eqLL$lon
    
EQ = AQ$EQ
   
###       /* calculate preliminary free solution */

  ###  maxITER = 5, residuals are eliminated between iterations
    AQ = XYlocate(Ldat,EQ,vel, 
      maxITER = maxITER[2],
      distwt = distwt,
      lambdareg =lambdareg ,
      FIXZ =FALSE,
      REG = REG,
      WTS = WTS,
      STOPPING = STOPPING,
      RESMAX =RESMAX ,
      tolx =   tolx,
      toly = toly ,
      tolz = tolz, PLOT=PLOT)


     Ksolutions[[2]] = AQ$guesses

       if(is.null(AQ$EQ))
       {
         return(list(EQ=NULL ))
         
         
       }
       
       

      eqLL = GEOmap::XY.GLOB(AQ$EQ$x, AQ$EQ$y, proj)
    AQ$EQ$lat = eqLL$lat
    AQ$EQ$lon = eqLL$lon
    
    EQ = AQ$EQ
    
    
###     /* calculate final free solution */
 ###   maxITER = 7
    AQ = XYlocate(Ldat,EQ,vel, 
      maxITER = maxITER[3],
      distwt = distwt,
      lambdareg =lambdareg ,
      FIXZ = FALSE,
      REG = REG,
      WTS = WTS,
      STOPPING = STOPPING,
       RESMAX = c(0,0),
      tolx =   tolx,
      toly = toly ,
      tolz = tolz, PLOT=PLOT)

   Ksolutions[[3]] = AQ$guesses

     if(is.null(AQ$EQ))
       {
         return(list(EQ=NULL ))
 

       }
       

###   /* if solution bad, try recovery */
   eqLL = GEOmap::XY.GLOB(AQ$EQ$x, AQ$EQ$y, proj)
    AQ$EQ$lat = eqLL$lat
    AQ$EQ$lon = eqLL$lon
    
    EQ = AQ$EQ
  
    KITS = AQ$its


    
    if(BADSOLUTION)
      {
       ###  maxITER = 4
        
        AQ = XYlocate(Ldat,EQ,vel, 
          maxITER = maxITER[4],
          distwt = distwt,
          lambdareg =2*lambdareg ,
          FIXZ=FALSE,
          REG = REG,
          WTS = WTS,
          STOPPING = STOPPING,
          tolx =   tolx,
          toly = toly ,
          tolz = tolz, PLOT=PLOT)


        eqLL = GEOmap::XY.GLOB(AQ$EQ$x, AQ$EQ$y, proj)
        AQ$EQ$lat = eqLL$lat
        AQ$EQ$lon = eqLL$lon
         EQ = AQ$EQ
         Ksolutions[[4]] = AQ$guesses

      }

########  finally wrap up and get error bars
    wup = eqwrapup(Ldat, EQ, vel, distwt=distwt, verbose=FALSE)
    
    return(list(EQ=EQ, ERR=wup, its=AQ$its ,  proj=proj, Ksolutions=Ksolutions ))
    

  }
