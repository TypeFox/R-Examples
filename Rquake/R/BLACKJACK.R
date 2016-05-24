BLACKJACK<-function(Ldat, vel)
  {
    
    
##########   jackknife an earthquake location
    usta = unique(Ldat$name)
    distwt = MeanStaDist(Ldat)

#####################   get location with all stations
    wstart = which.min(Ldat$sec)

#########   choose initial guess of earthquake beneath earliest arrival stations
    EQ = list(lat=Ldat$lat[wstart], lon=Ldat$lon[wstart], z=6, t=Ldat$sec[wstart] )
    
##########   locate earthquake
    AQ = Vlocate(Ldat,EQ,vel,
      distwt = distwt,
      lambdareg =100 ,
      REG = TRUE,
      WTS = TRUE,
      STOPPING = TRUE,
      tolx =   0.01,
      toly = 0.01 ,
      tolz = 0.05,
      maxITER = c(17,5,17,4) ,
      RESMAX = c(0.1, 0.1),
      PLOT=FALSE)
    
    LAM = c(AQ$EQ$lat, AQ$EQ$lon, AQ$EQ$z)
#########################

#########################
#########################
#########################    get psuedo values
    JQ = vector(mode="list")
    
    for(i in 1:length(usta))
      {
        LeaveOut = usta[i]

        wout =  which( Ldat$name!=LeaveOut )
        
        Ddat =  data.frame(Ldat)
        Ddat =  Ddat[wout , ]
        Adat = as.list(Ddat)

        wstart = which.min(Adat$sec)

#########   choose initial guess of earthquake beneath earliest arrival stations
        EQ = list(lat=Adat$lat[wstart], lon=Adat$lon[wstart], z=6, t=Adat$sec[wstart] )

        JQ[[i]] = Vlocate(Adat,EQ,vel,
            distwt = distwt,
            lambdareg =100 ,
            REG = TRUE,
            WTS = TRUE,
            STOPPING = TRUE,
            tolx =   0.01,
            toly = 0.01 ,
            tolz = 0.05,
            maxITER = c(17,5,17,4) ,
            RESMAX = c(0.1, 0.1),
            PLOT=FALSE)


      }

    BEDS = matrix(ncol=3, nrow=length(JQ))
    for(i in 1:length(usta))
      {
        BEDS[i, ] = c(JQ[[i]]$EQ$lat,JQ[[i]]$EQ$lon, JQ[[i]]$EQ$z )

      }
    rownames(BEDS)<-usta

  ##  plot(GX, GY, type='n', xlab="km" , ylab="km" , asp=1)
  ##  points(XY, pch=6)


    BEXY = GEOmap::GLOB.XY(BEDS[,1], BEDS[,2], AQ$proj)
   ## points(BEXY, pch=8, col='red')


    sweep(BEDS,2, LAM)

    ##   bi = (n-1)*(lam- lami)
    ##  or  bi = (1-n)*(lami-lam)

    #######  station influence
    SI = (1-length(usta))*sweep(BEDS,2, LAM)

   #######   BIAS = -apply(BIpsuedo, 2, sum)/length(usta)

    LSI = list(u=usta,SI=SI,BEXY=BEXY , BEDS=BEDS, AQ=AQ  ) 
    return( LSI )
    
  }
##################
##################
##################


