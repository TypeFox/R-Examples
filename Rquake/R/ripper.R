ripper<-function(AQ)
  {

    zip = AQ$Ksolutions

    iz = 1
    miz =  zip[[iz]]
    co = rep(iz, length(miz[,1]))
    
    rat = cbind(  cbind(miz, co))
    
    
    for(iz in 2:length(zip))
      {
        miz =  zip[[iz]]
        co = rep(iz, length(miz[,1]))
        
        rat = rbind(rat ,  cbind(miz, co))
        
      }

    gat = GEOmap::XY.GLOB(rat[,1], rat[,2], AQ$proj)
    
    rat[,1:2] = cbind(gat$lat, gat$lon)
    return(rat)
    
  }



