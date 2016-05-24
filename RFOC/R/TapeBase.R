TapeBase<-function( )
  {
###   TATE plot from Tape and Tape 

  ###  require(GEOmap)
   
    d1=list()
    d1$'lam'=c(0,0,0,-30,-30,-30,-30,30,30,30,30)
    d1$'phi'=c(-90,0,90,60.5038,35.2644,0,-54.7356,-60.5038,-35.2644,0,54.7356)
    d1$'name'=c("ISO", "DC", "ISO", "C", "LVD", "CLVD", "-", "C", "LVD", "CLVD", "-")
    d1$'a1'=c(0,15,0,-18,-22.5,-30,-15,18,22.5,30,15)
    d1$'a2'=c(-15,15,15,0,0,0,0,0,0,0,0)

###########  first PATCH
    g = GEOmap::getgreatarc(d1$phi[5],d1$lam[5],d1$phi[11],d1$lam[11],    100)
    
    PATCH11 = HAMMERprojXY(g$lat*pi/180, g$lon*pi/180)

    lats = seq(from=90, to=d1$phi[5], length=20)
    lons = rep(-30, length=20)
    
    PATCH12 = HAMMERprojXY(lats*pi/180, lons*pi/180)

    lats = seq(from=d1$phi[11], to=90, length=20)
    lons = rep(30, length=20)
    
    PATCH13 = HAMMERprojXY(lats*pi/180, lons*pi/180)

    POLYh1 = list(x=NULL, y=NULL)
    POLYh1$x =c(PATCH12$x, PATCH11$x, PATCH13$x )
    
    POLYh1$y =c(PATCH12$y, PATCH11$y, PATCH13$y)



###########  second PATCH
    g = GEOmap::getgreatarc(d1$phi[7],d1$lam[7],d1$phi[9],d1$lam[9],    100)
    
    PATCH21 = HAMMERprojXY(g$lat*pi/180, g$lon*pi/180)

    lats = seq(to=-90, from=d1$phi[9], length=20)
    lons = rep(30, length=20)
    
    PATCH22 = HAMMERprojXY(lats*pi/180, lons*pi/180)

    lats = seq(from=-90, to=d1$phi[7], length=20)
    lons = rep(-30, length=20)
    
    PATCH23 = HAMMERprojXY(lats*pi/180, lons*pi/180)

    POLYh2 = list(x=NULL, y=NULL)
    POLYh2$x =c(PATCH21$x, PATCH22$x, PATCH23$x )
    
    POLYh2$y =c(PATCH21$y, PATCH22$y, PATCH23$y)

    ##   polygon(h2)

    
    lons = seq(from=-30, to=30, by=10)*pi/180
    
    lats = seq(from=-90, 90, by=10)*pi/180

    
    Left1 =  HAMMERprojXY(0, min(lons))
    right2 =  HAMMERprojXY(0, max(lons))
    
    top1 =  HAMMERprojXY(lats[length(lats)], 0)
    bot1 =  HAMMERprojXY(lats[1], 0)

    ## Underneath: filled Polygon patches 
   ## polygon(POLYh1$x, POLYh1$y, col=pcol[1] , border=NA)
    
  ##  polygon(POLYh2$x, POLYh2$y, col=pcol[2] , border=NA)
    

    ####  longitude and latitude lines
    LONSp1 = vector(mode="list")
 LONSp2  = vector(mode="list")
LATSp1  = vector(mode="list")
      for(i in 1:length(lats))
      {
        
        LONSp1[[i]] =  HAMMERprojXY(lats[i], min(lons))
        LONSp2[[i]] =  HAMMERprojXY(lats[i], max(lons))
        
      }
    
    for(i in 1:length(lons) )
      {
        
        LATSp1[[i]] =  HAMMERprojXY(lats, lons[i] )
       
        
      }

    ###

    PTSh1 = HAMMERprojXY(d1$phi*pi/180, d1$lam*pi/180)

#####  horizontal zero line
    HOZk1=list()
    HOZk1$'lam'=seq(from=-30, to =30, length=20)
    HOZk1$'phi' = rep(0, length=length(HOZk1$'lam'))

    
    HOZh1 = HAMMERprojXY(HOZk1$phi*pi/180, HOZk1$lam*pi/180)
   
###  vertical zero line
    VERTk1=list()
    VERTk1$'phi'=seq(from=-90, to=90, length=2)
    VERTk1$'lam'= rep(0, length=length(VERTk1$'phi'))

 VERTh1 = HAMMERprojXY(VERTk1$phi*pi/180, VERTk1$lam*pi/180)
   
 ###        lines(VERTh1$x, VERTh1$y, lty=2, lwd=2)

###  add in a bold dotted line for LVD-1
 ###    lines(PATCH11$x, PATCH11$y, lty=2, lwd=2)
###  add in a bold dotted line for LVD-2
    
 ###    lines(PATCH21$x, PATCH21$y, lty=2, lwd=2)
###  connect up the Crack lines
    CRACKh1 = HAMMERprojXY(d1$phi[c(4,8) ]*pi/180, d1$lam[c(4,8) ]*pi/180)
###     lines(CRACKh1$x, CRACKh1$y, lty=2, lwd=2)

BLIST = list(d1, Left1, right2, top1, bot1, POLYh1, POLYh2,LONSp1, LONSp2,LATSp1,
  PTSh1, HOZh1, VERTh1, PATCH11, PATCH21, CRACKh1)
names(BLIST) <-c("d1", "Left1", "right2", "top1", "bot1", "POLYh1", "POLYh2",
                 "LONSp1", "LONSp2","LATSp1",
  "PTSh1", "HOZh1", "VERTh1", "PATCH11", "PATCH21", "CRACKh1")
    
  invisible(BLIST)
    

  }
