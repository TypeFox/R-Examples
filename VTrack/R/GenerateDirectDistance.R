GenerateDirectDistance <-
function(sPointsFile)
{
    Latitude <- sPointsFile[,2]	
    Longitude <- sPointsFile[,3]
    LatLon <- cbind(Latitude, Longitude)*pi/180	
    rDetectionRadius <- sPointsFile[,4]/1000
    nxy <- length(unique(sPointsFile[,1]))
    d <- matrix(0,nrow=nxy, ncol=nxy)
    dr <- matrix(0,nrow=nxy, ncol=nxy)
    
    # generates distance matrix
    for(i in 1:(nxy-1)) {		
      d[i,i] <- 1.0
      for (j in (i+1):nxy) {			
        d[i,j] <- sin(LatLon[i,1])*sin(LatLon[j,1]) + cos(LatLon[i,1]) * cos(LatLon[j,1])*cos(abs(LatLon[i,2] - LatLon[j,2]))			
        d[j,i] <- d[i,j]		
      }	
    }
    
    d[nxy, nxy] <- 1.0
    d <- acos(d)
       
    for (i in 1: nxy) { 	
      for (j in i:nxy) { 		
        if (d[i,j] < 0.00000000001) {d[i,j] <- 0}
        d[j,i] <- d[i,j]	 	
      }
    }
    d <- ifelse (d < 0.000000000001, 0.0, d)
    d <- 6371*d
    
    # generates detection radius matrix
    for(j in 1:nxy)
    {
      for (i in j:nxy)
      {
        if(i != j)
        {
          #summate the distance between i & j
          dr[i,j] <- rDetectionRadius[j] + rDetectionRadius[i]
          dr[j,i] <- dr[i,j]
        }
      }
    }
    # subtract detection radius from overall distance
    d1 <- d-dr
    # Make sure no negative distances
    for(i in 1:length(d1)) if(d1[i]<=0) d1[i]=0 else  d1[i]=d1[i]
    d1 <- data.frame(DM=sPointsFile[,1],d1)
    # Assign receiver/station names  
    colnames(d1) <- c("DM",sPointsFile[,1])
  return(d1)
}
