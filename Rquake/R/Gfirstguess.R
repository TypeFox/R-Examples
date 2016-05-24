Gfirstguess<-function(Ldat, type="first" )
  {

    if(missing(type)) { type="first" }
    
    if(identical(type, "first"))
      {
        A1T = Qrangedatetime(Ldat)
        wmin = which.min(Ldat$sec)
        pop = c(Ldat$lat[wmin],Ldat$lon[wmin],Ldat$z[wmin],Ldat$sec[wmin] )
      }
    
    if(identical(type, "mean"))
      {
        
        MLAT = mean(Ldat$lat)
        MLON = mean(Ldat$lon)
        MZ = mean(Ldat$z)
        wmin = 1
        pop = c(MLAT,MLON, MZ,Ldat$sec[wmin] )
        
      }
  if(identical(type, "median"))
      {
        
        MLAT = median(Ldat$lat)
        MLON = median(Ldat$lon)
        MZ = median(Ldat$z)
        wmin = 1
        pop = c(MLAT,MLON, MZ,Ldat$sec[wmin] )
        
      }
    
    return(pop) 
  }
