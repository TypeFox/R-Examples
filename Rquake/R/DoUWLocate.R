DoUWLocate<-function(LF, stas, vel,
                    params=list(distwt = 100,
          lambdareg=20,
          REG = TRUE,
          WTS = TRUE,
          STOPPING = TRUE,
          tolx = 0.005,
          toly = 0.005,
          tolz = 0.01))
  {

    Gback = vector(mode="list")
    
    for(i in 1:length(LF))
      {
        g1 = RSEIS::getpfile( LF[i] )

        m1 = match(g1$STAS$name, stas$name)

        g1$STAS$lat = stas$lat[m1]
        g1$STAS$lon = stas$lon[m1]
        g1$STAS$z = stas$z[m1]


       
        ##  points(g1$H$lon, g1$H$lat, pch=8, col='red')

        w1 = which(!is.na(g1$STAS$lat))
        sec = g1$STAS$sec[w1]

        N = length(sec)


        

####   set up the data list for earthquake location
        
        Ldat =    list(
          name = g1$STAS$name[w1],
          sec = g1$STAS$sec[w1],
          phase = g1$STAS$phase[w1],
          lat=g1$STAS$lat[w1],
          lon = g1$STAS$lon[w1],
          z = g1$STAS$z[w1],
          err= g1$STAS$err[w1],
          yr = rep(g1$LOC$yr , times=N),
          jd = rep(g1$LOC$jd, times=N),
          mo = rep(g1$LOC$mo, times=N),
          dom = rep(g1$LOC$dom, times=N),
          hr =rep( g1$LOC$hr, times=N),
          mi = rep(g1$LOC$mi, times=N) )





        Ldat$err[Ldat$err<=0] = 0.05
        
        Ksta = length( unique(Ldat$name) )
        
        
        
###  match up picks with stations from the station file
        
        ####   set up the data list for earthquake location
        cat(paste("#################      ", i, Ksta), sep="\n" )
#### print( data.frame(Ldat) )
        Ldat = LeftjustTime(Ldat)
        
        MinDate = list(yr=Ldat$yr[1], jd=Ldat$jd[1],  mo=Ldat$mo[1],  dom=Ldat$dom[1], 
          hr=Ldat$hr[1], mi=Ldat$mi[1], sec=0)
        
        MLAT = median(Ldat$lat)
        MLON = median(Ldat$lon)

        proj = GEOmap::setPROJ(type=2, LAT0=MLAT, LON0=MLON)

        XY = GEOmap::GLOB.XY(Ldat$lat, Ldat$lon, proj)

        wstart = which.min(Ldat$sec)

        EQ = list(lat=Ldat$lat[wstart], lon=Ldat$lon[wstart], z=6, t=0)

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
          tolz =  params$tolz
          )

        MinDate$sec = zypo$EQ$t
        zypo$EQ$Time = MinDate
        Gback[[i]] = zypo

      }

    return(Gback)
  }
