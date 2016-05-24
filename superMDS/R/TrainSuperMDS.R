TrainSuperMDS <-
function(d=NULL, y, alpha=.5, S=2, x=NULL, nstarts=5, silent=FALSE){
  outputs <- list()
  crits <- rep(NA, nstarts)
  if(!silent) cat("Starting iter 1", fill=TRUE)
  outputs[[1]] <- TrainSuperMDSOnce(d=d,y=y,alpha=alpha,S=S,x=x,z=NULL)
  crits[1] <- min(outputs[[1]]$crits)
  if(nstarts>1){
    for(iter in 2:nstarts){
      if(!silent) cat("Starting iter", iter, fill=TRUE)
      z <- matrix(rnorm(nrow(outputs[[1]]$d)*S), nrow=nrow(outputs[[1]]$d))
      z[y==2,] <- z[y==2,]+5
      outputs[[iter]] <- TrainSuperMDSOnce(d=d,y=y,alpha=alpha,S=S,x=x,z=z)
      crits[iter] <- min(outputs[[iter]]$crits)
    }
  }
  return(outputs[[which.min(crits)]])
}


 
