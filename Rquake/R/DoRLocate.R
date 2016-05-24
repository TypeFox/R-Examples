DoRLocate<-function(LF, stas, vel,
                    params=list(distwt = 100,
                      lambdareg=20,
                      REG = TRUE,
                      WTS = TRUE,
                      STOPPING = TRUE,
                      tolx = 0.005,
                      toly = 0.005,
                      tolz = 0.01,
                      RESMAX = c(4,5),
                      maxITER = c(7, 5, 7, 4)
                      ))
  {

    Gback = vector(mode="list")
    twpx = list()
    
    for(i in 1:length(LF))
      {

        load(LF[i])
        
        if(any(twpx$onoff<1))
          {
            Ldat = RSEIS::deleteWPX(twpx, which(twpx$onoff<1) )
          }
        else
          {
            Ldat = twpx 
          }

        if(is.null(Ldat)) next
        
        Ldat$err[Ldat$err<=0] = 0.05

        Ksta = length( unique(Ldat$name) )

        
        maxITER = params$maxITER
        if(Ksta<4) {  maxITER = c(7, 0,0 , 4) }
        ###  match up picks with stations from the station file

        ####   set up the data list for earthquake location
        cat(paste("#################      ", i, Ksta), sep="\n" )
        #### print( data.frame(Ldat) )
        Ldat = LeftjustTime(Ldat)

        MinDate = list(yr=Ldat$yr[1], jd=Ldat$jd[1],  mo=Ldat$mo[1],  dom=Ldat$dom[1], 
                        hr=Ldat$hr[1], mi=Ldat$mi[1], sec=0)
          
     m1 = match(Ldat$name, stas$name)

        Ldat$lat = stas$lat[m1]
        Ldat$lon = stas$lon[m1]
        Ldat$z = stas$z[m1]

       
        MLAT = median(Ldat$lat)
        MLON = median(Ldat$lon)

        proj = GEOmap::setPROJ(type=2, LAT0=MLAT, LON0=MLON)

        XY = GEOmap::GLOB.XY(Ldat$lat, Ldat$lon, proj)

        wstart = which.min(Ldat$sec)

        EQ = list(lat=Ldat$lat[wstart], lon=Ldat$lon[wstart], z=6, t=Ldat$sec[wstart]-0.1 )

        eqxy = GEOmap::GLOB.XY(EQ$lat, EQ$lon, proj)

        EQ$x = XY$x[wstart]
        EQ$y = XY$y[wstart]

        zypo = Vlocate(
          Ldat,EQ,vel, 
          distwt = params$distwt,
          lambdareg= params$lambdareg,
          REG =  params$REG,
          WTS =  params$WTS,
          STOPPING = params$STOPPING,
          tolx =  params$tolx,
          toly =  params$toly,
          tolz =  params$tolz,
          maxITER = maxITER
          )

        if(is.null(zypo$EQ))
          {
            
            zypo = list(EQ=list(lat=NA, lon=NA), Time=MinDate)
          }

        MinDate$sec = zypo$EQ$t
        zypo$EQ$Time = MinDate

        zypo$vel = vel
        zypo$PickTimes = Ldat
        Gback[[i]] = zypo


        
      }

    return(Gback)
  }
