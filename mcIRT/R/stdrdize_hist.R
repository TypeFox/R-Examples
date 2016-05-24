stdrdize_hist <- function(At,quadnodes)
{
  
  weigm <- weighted.mean(quadnodes,At) # compute weighted mean
  quadnodesC <- quadnodes - weigm      # centering nodes
  
  weigSD <- sqrt(sum((At*(quadnodesC)^2)/sum(At))) # compute weighted sd
  
  quadnodesCS <- quadnodesC/weigSD # scale nodes
  
  At_new <- At *weigSD  # rescale weights
  
  DELTA <- quadnodesCS[2] - quadnodesCS[1] # estimate the difference between two adjacent categories
  LG <- length(At_new)  # number of nodes
  
  expolERG <- sapply(quadnodes,function(x) # extra/interpolate between nodes - see woods article of EH
  {
    if(x < min(quadnodesCS))
    {
      
      eins <- (At_new[1]/At_new[2])^((quadnodesCS[1] - x)/DELTA)
      Nq <- eins * At_new[1]
      
    } else if(x > max(quadnodesCS))
    {
      eins <- (At_new[LG]/At_new[LG-1])^((x - quadnodesCS[LG])/DELTA)
      Nq <- eins * At_new[LG]
    } else {
      
      drueber <- which(x < quadnodesCS)[1] # gibt mir den darüberliegenden node
      drunter <- which(x > quadnodesCS)[length(which(x > quadnodesCS))] # gibt mir jenen node der strikt kleiner ist als mein dings     
      Nq <- (x - quadnodesCS[drunter])/DELTA * (At_new[drueber] - At_new[drunter]) + At_new[drunter]
    }
    
  }) 

  # normalize and name
  weights <- expolERG / sum(expolERG)
  nodes   <- quadnodes
  
  return(list(nodes=nodes,weights=weights))  
   
}