DistIdealPatt <-
function(Y,Q,weight){
  
   
  #-----------------------Basic variables----------------------#
  #N:number of examinees
  #J:number of items
  #K:number of attributes
  #M:number of ideal attribute patterns, which is equal to 2^K
  N <- dim(Y)[1]
  J <- dim(Y)[2]
  K <- dim(Q)[2]
  M <- 2^K
  #------------------------------------------------------------#
  #A:alpha matrix
  #E:eta matrix
  
  A <- alpha (K)
  E <- eta (K,J,Q)
  
  #=============================Compute weight===================================#
  if (is.null(weight)){ #unspecified weights
    p <- apply(Y,2,mean)
    
    if (min(p) == 0 | max(p) == 1)
    {
      warning("Cannot compute weights because some weights equal to NA or Inf, unweighted Hamming distance will be used.")
      weight <- c(rep(1,times = J))
    }else{
      weight <- 1/(p*(1-p))
    }
  }
  #============================compute distance===================================#
  dis <- c(rep(NA,M))
  for (i in 1:M){
    dis[i] <- apply(as.matrix(apply(abs(Y-E[rep(i,times = N),]),2,sum)*weight),2,sum)/N
  }
  output <- list(dist = dis,weight = weight)
  return(output)
}
