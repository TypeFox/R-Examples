RElocateEQ<-function(lps, sta, vel, cont=TRUE,
                     sleep=0.5, mapfun=NULL, PLOT=TRUE, proj=NULL, xlim=NULL, ylim=NULL)
  {
##################
############    given a set of pickfiles,
###########   show relocations

    if(missing(mapfun)) { mapfun=NULL  }

  #####  whatmap = mapfun

    
    if(!is.null(mapfun))
      {
   #####    ### print("assigning whatmap" )
        assign("whatmap", mapfun )
        
      }
    else
     {
        whatmap = NULL

      }
    
    grcol = grey(seq(from=0.3, to=0.95, length=50))
    
    Ntot = length(lps)

    Qout  = vector(mode="list")
    

    
    for(ipi in 1:Ntot)
      {
        g1 = RSEIS::getpfile(lps[ipi], stafile = NULL)

        MA = match(g1$STAS$name, sta$name)

        g1$STAS$lat = sta$lat[MA]
        g1$STAS$lon = sta$lon[MA]
        g1$STAS$z = sta$z[MA]

        w1 = which(!is.na(g1$STAS$lat) & !is.na(g1$STAS$lon) )
        Ldat  = LDATlist(g1, w1)

        
        distwt = MeanStaDist(Ldat)

        
#############   get earliest arrival time
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
          tolz = 0.05, maxITER = c(7,5,7,4) , RESMAX = c(0.1, 0.1),  PLOT=FALSE)

        
        if(PLOT==TRUE)
          {
            plotEQ(Ldat, AQ, add=FALSE, prep=TRUE, proj=proj, xlim=xlim, ylim=ylim )


            if(cont)
              {
                contPFarrivals(g1, sta, proj=AQ$proj,cont=TRUE, POINTS=FALSE,
                               image=TRUE , col=grcol,     phase="P", add=TRUE)
              }

            if(!is.null(whatmap))
              {
### print("mapping")
                whatmap(AQ$proj)

              }

            
            plotEQ(Ldat, AQ, add=TRUE)
            Sys.sleep(sleep)
          }
        Qout[[ipi]] = AQ
        
      }

    invisible(Qout)
  }

