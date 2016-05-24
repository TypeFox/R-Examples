
.runLengthShroNorm <- function(x, mean, sigma, n, delta, ubd){
  
  limit <- length(x)
  
  m = 2 
  wm = 0
  wmv<- rep(NA, limit)
  
  while(m < limit && wm < ubd){
    s1 = 0
    wm = 0 
    
    for(i in 1:(m - 1)){
      s1 = s1 + x[m - i + 1] - mean
      wm = wm + exp(-i * n * (delta ^ 2) / (2 * sigma ^ 2) + n * delta * s1 / sigma ^ 2)
    }
    
    wmv[m] <- wm
    
    if(wm > ubd || (m+1) == limit){
      res <- vector("list", 0)
      if(wm > ubd){
        res$rl <- m
        res$w <- wmv[1:m]
      } else {
        res$rl <- Inf
        res$w <- wmv
      }
    }   
    m = m + 1
  }
  return(res)
}


.runLengthShroPois <- function(x, rho, delta, ubd){
  
  limit <- length(x)
  
  m = 2 
  wm = 0
  wmv<- rep(NA, limit)
  
  while(m < limit && wm < ubd){
    s1 = 0
    wm = 0 
    
    for(i in 1:(m - 1)){
      s1 = s1 + x[m - i + 1]
      wm = wm + exp(-i * delta + s1 * log(rho))
    }
    
    wmv[m] <- wm
    
    if(wm > ubd || (m+1) == limit){
      res <- vector("list", 0)
      if(wm > ubd){
        res$rl <- m
        res$w <- wmv[1:m]
      } else {
        res$rl <- Inf
        res$w <- wmv
      }
    }   
    m = m + 1
  }
  return(res)
}


shroArlPfaCedNorm <- function (mean0=0, mean1=NA, sd=1,
                     n = 10, delta=1,
                     tau=NA,
                     N=100, limit=10000, seed=NA, 
                     w=19, 
                     printSummary=TRUE) 
{
  if(is.na(delta))
    delta = (mean1 - mean0)/sd
  
  if(is.na(mean1))
    mean1 = mean0 + delta*sd
  
  if(!is.na(seed))
    set.seed(seed)
  
  if(is.na(tau)){
    data <- matrix(rnorm(n=limit*N, mean0, sd/sqrt(n)), nrow=N)
  } else {
    data <- matrix(rnorm(n=(tau)*N, mean0, sd/sqrt(n)), nrow=N)
    data <- cbind(data, matrix(rnorm(n=(limit-tau)*N, mean1, sd/sqrt(n)), nrow=N))  
  }
  
  res <- list(run=apply(data, MARGIN=1, .runLengthShroNorm, mean=mean0, sigma=sd, n=n, delta=delta, ubd=w))
  res$rls <- sapply(res$run, function(x) x$rl)
  
  if(is.na(tau)){
    res$statistic <- c(mean(res$rls), 
             sqrt((mean(res$rls^2) - mean(res$rls)^2)/N))
    names(res$statistic) <- c("ARL", "Std. Error")
  } else {
    pfa <- sum(res$rls < tau)/nrow(data)
    ced = sum(res$rls[res$rls >= tau])/(nrow(data) - sum(res$rls < tau)) - tau
    se = sqrt((sum(res$rls[res$rls >= tau]^2)/(nrow(data) - sum(res$rls < tau)) - ced ^ 2) / (nrow(data) - sum(res$rls < tau)))
    res$statistic <- c(mean(res$rls), 
             sqrt((mean(res$rls^2) - mean(res$rls)^2)/N),
             pfa, ced, se)
    names(res$statistic) <- c("ARL", "Std. Error", "PFA", "CED", "CED-Std. Error")
  }
  
  if(printSummary)
    print(res$statistic)

  invisible(res)
}


shroArlPfaCedPois <- function (lambda0=10, lambda1=NA,
                               delta=1,
                               tau=NA,
                               N=100, limit=10000, seed=NA, 
                               w=19, 
                               printSummary=TRUE) 
{
  if(is.na(delta))
    delta <- (lambda1 - lambda0)
  
  if(is.na(lambda1))
    lambda1 <- lambda0 + delta
  
  rho <- (lambda0 + delta)/lambda0
  
  if(!is.na(seed))
    set.seed(seed)
  
  if(is.na(tau)){
    data <- matrix(rpois(n=limit*N, lambda0), nrow=N)
  } else {
    data <- matrix(rnorm(n=(tau)*N, lambda0), nrow=N)
    data <- cbind(data, matrix(rnorm(n=(limit-tau)*N, lambda1), nrow=N))  
  }
  
  res <- list(run=apply(data, MARGIN=1, .runLengthShroPois, rho=rho, delta=delta, ubd=w))
  res$rls <- sapply(res$run, function(x) x$rl)
  
  if(is.na(tau)){
    res$statistic <- c(mean(res$rls), 
             sqrt((mean(res$rls^2) - mean(res$rls)^2)/N))
    names(res$statistic) <- c("ARL", "Std. Error")
  } else {
    pfa <- sum(res$rls < tau)/nrow(data)
    ced = sum(res$rls[res$rls >= tau])/(nrow(data) - sum(res$rls < tau)) - tau
    se = sqrt((sum(res$rls[res$rls >= tau]^2)/(nrow(data) - sum(res$rls < tau)) - ced ^ 2) / (nrow(data) - sum(res$rls < tau)))
    res$statistic <- c(mean(res$rls), 
             sqrt((mean(res$rls^2) - mean(res$rls)^2)/N),
             pfa, ced, se)
    names(res$statistic) <- c("ARL", "Std. Error", "PFA", "CED", "CED-Std. Error")
  }
  
  if(printSummary)
    print(res$statistic)
  
  invisible(res)
}


