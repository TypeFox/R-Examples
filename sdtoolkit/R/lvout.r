#a helper function for null probability - finds subsets of datasets based on restrictions

lvout <- function(dset,y=NULL,lbox){

  rmdsizes <- vector(length = length(lbox[[1]]))
  rmdhighs <- vector(length = length(lbox[[1]]))

  for (i in 1:length(lbox[[1]])){
    tempbox <- list()
    
    tempbox[[1]] <- lbox[[1]][-i]
    tempbox[[2]] <- lbox[[2]][-i,]     
    
    if (is.null(y)){  
      allin  <-  unionpts(dset,list(tempbox)) 
    } else{
      allin  <-  unionpts(cbind(dset,y),list(tempbox))
    }
      
    rmdsizes[i]  <- sum(allin)
    if(is.null(y)){
      rmdhighs[i]  <- sum(dset[allin,ncol(dset)])
    } else{
      rmdhighs[i]  <- sum(y[allin])
    }
  }
  
  rstats <- cbind(rmdsizes,rmdhighs,rmdhighs/rmdsizes)

  #While at it, get the original stats and return that too
  tempbox <- lbox

    if (is.null(y)){  
      allin  <-  unionpts(dset,list(tempbox)) 
    } else{
      allin  <-  unionpts(cbind(dset,y),list(tempbox))
    }

    origtotin <- sum(allin)
  
    if(is.null(y)){
      orighighin  <- sum(dset[allin,ncol(dset)])
    } else{
      orighighin  <- sum(y[allin])
    }

  attr(rstats,"origtotin") <- origtotin
  attr(rstats,"orighighin") <- orighighin
  attr(rstats,"origdens") <- orighighin/origtotin

  return(rstats)

}

