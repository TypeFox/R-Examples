w.matrix <- function(Coords, Bandwidth, WType='Binary', family='adaptive'){
  
    Distances <- dist(Coords)
    Dij <- as.matrix(Distances)
    
    Obs <- length(Coords[,1])
    
    if(family =='adaptive' && Bandwidth>=Obs){
      Bandwidth<-Obs-1
      msg<-cat("The number of nearest neighbours exceeds the number of observations minus one.\nBandwidth set to:", Bandwidth)
    }
    
    Wts <- matrix(data=0, nrow=Obs, ncol=Obs)
    
    for(i in 1:Obs)
      {
        #Create a dataset of the distances 
        DNeighbour <- Dij[,i]
      
        #Sort by distance
        DSorted <- sort(DNeighbour)
      
        if(family=='adaptive')
          { 
            #Keep Nearest Neighbours
            SubSet1 <- DSorted[2:(Bandwidth+1)]
          
            #Find furthest neighbour
            Kernel_H <- max(SubSet1)
          } 
        else 
          { 
            if(family=='fixed')
              {
                Kernel_H <- Bandwidth
              }
          }
        
      #Calculate weights
      
      for(j in 1:Obs)
        {
          if (DNeighbour[[j]] > Kernel_H)
            {
              Wts[i,j]<-0 
            }
          else
            {
              if(WType=='Bi-square')
                {
                  Wts[i,j]<-(1-(DNeighbour[[j]]/Kernel_H)^2)^2
                }
              else
                {
                  Wts[i,j]<-1
                }
          
            }
         if (j==i)
           {
              Wts[i,j]<-0
           }
        
        }
      
        if(WType=='RSBi-square')
          {
            Wts[i,] <- Wts[i,]/sum(Wts[i,])
          }
      
      }
    
      return(Wts)
}
    
  