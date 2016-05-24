EvalOneBlock <- function(da, alpha, tol, nbinit){
  if (length(alpha)>1){
    tmp <- uniquecombs(da)
    weight <- as.numeric(table(attr(tmp,"index")))
    bic <- XEMblock(as.matrix(tmp), weight, alpha, tol, nbinit)$loglike - length(alpha) * log(sum(weight))
  }else{
    bic <- log(alpha)*sum(da) + log(1-alpha)*sum(1-da) - 0.5*log(sum(da) + sum(1-da))
  }
  return(bic)
}

OneMH <- function(x, alpha, tol, nbinit, iterMH){
  # initalisation
  d <- ncol(x)
  omega_best <- omega_cand <- omega_current <- sample(1:d,d,replace=TRUE)
  values_best <-  rep(0, d)
  for (b in unique(omega_best))  values_best[b] <- EvalOneBlock(x[,which(omega_best==b)], alpha[which(omega_best==b)], tol, nbinit)
  values_cand <- values_current <- values_best
  
  iter <- 0
  
  while (iter<iterMH){
    iter <- iter + 1
    # generation du candidat
    omega_cand <- omega_current
    jmove <- sample(1:d,1)
    blockdepart <- omega_current[jmove]
    if (any(omega_current==0)){
      blockarrive <- sample(c(unique(omega_current),which(values_current==0)[1]),1)
    }else{
      blockarrive <- sample(1:d,1)
    }
    omega_cand[jmove] <- blockarrive
    values_cand <- values_current
    if (any(omega_cand==blockdepart)){
      values_cand[blockdepart] <- EvalOneBlock(x[,which(omega_cand==blockdepart)], alpha[which(omega_cand==blockdepart)], tol, nbinit)
    }else{
      values_cand[blockdepart] <- 0
    }    
    values_cand[blockarrive] <- EvalOneBlock(x[,which(omega_cand==blockarrive)], alpha[which(omega_cand==blockarrive)], tol, nbinit)
    rho <- exp(sum(values_cand[blockarrive] + values_cand[blockdepart] - values_current[blockarrive] - values_current[blockdepart]))  
    if( runif(1)<rho){
      values_current <- values_cand
      omega_current <- omega_cand
      if (sum(values_current)>sum(values_best)){
        values_best <- values_current
        omega_best <- omega_current
        iter <- 0
      }
    }
  }
  outomega <- rep(0, length(omega_best))
  uni <- unique(omega_best)
  for (v in 1:length(uni)) outomega[which(omega_best==uni[v])] <- v
  return(list(blocks=outomega, bic=sum(values_best)))  
}


