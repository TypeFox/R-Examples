
.runLength <- function(x, kp, km, ubd, lbd, side){
  
  limit <- length(x)
  undercontrol <- NA
  cusum.pos <- rep(NA, limit)
  if(side %in% c("both", "upper")){
    z.p <- x - kp
    cusum.pos[1] <- max(0, z.p[1])
    if(cusum.pos[1] > ubd)
      undercontrol <- FALSE
    else
      undercontrol <- TRUE
    i <- 2
    
    while(undercontrol && i <= limit){
      cusum.pos[i] <- max(0, cusum.pos[i - 1] + 
                            z.p[i])
      if(cusum.pos[i] > ubd)
        undercontrol <- FALSE
      i <- i +1
    }
  }
  
  cusum.neg <- rep(NA, limit)
  if(side %in% c("both", "lower")){
    z.m <- x - km
    cusum.neg[1] <- min(0, z.m[1])
    if(is.na(undercontrol) || i < limit)
      limit <- i
    if(cusum.neg[1] < lbd)
      undercontrol <- FALSE
    else
      undercontrol <- TRUE
    i <- 2
    
    while(undercontrol && i <= limit){
      cusum.neg[i] <- min(0, cusum.neg[i - 1] + 
                            z.m[i])
      if(cusum.neg[i] < lbd)
        undercontrol <- FALSE
      i <- i +1
    }
  }
  
  violations <- list(lower = which(cusum.neg < lbd), 
                     upper = which(cusum.pos > ubd))
  res <- vector("list", 0)
  res$violationsLower <- which(cusum.neg < lbd)
  res$violationsUpper <- which(cusum.pos > ubd)
  if(length(res$violationsLower) == 0)
    res$violationsLower <- Inf
  if(length(res$violationsUpper) == 0)
    res$violationsUpper <- Inf
  
  res$rl <- switch(side,
                   both = min(c(res$violationsLower, res$violationsUpper)),
                   upper= min(res$violationsUpper),
                   lower=min(res$violationsLower))
  return(res)
}

.cusumPfaCed <- function(data, tau, kp=1, km=-1, hp=3, hm=-3, side="both", printSummary){
  
  if (!side %in% c("both", "upper", "lower")) 
    warning(gettextf("side = '%s' is not supported. Using 'both'", 
                     side), domain = NA)
  
  n <- ncol(data)
  
  res <- list(run = apply(data, MARGIN=1, .runLength, kp=kp, km=km, ubd=hp, lbd=hm, side=side))
  res$rls <- sapply(res$run, function(x) x$rl)
  
  pfa <- sum(res$rls < tau)/nrow(data)
  ced = sum(res$rls[res$rls >= tau])/(nrow(data) - sum(res$rls < tau)) - tau
  se = sqrt((sum(res$rls[res$rls >= tau]^2)/(nrow(data) - sum(res$rls < tau)) - ced ^ 2) / (nrow(data) - sum(res$rls < tau)))
  
  res$statistic <- c(pfa, ced, se)
  names(res$statistic) <- c("PFA", "CED", "Std. Error")
  
  if(printSummary)
    print(res$statistic)

  invisible(res)
}

cusumArl <- function (..., randFunc = rnorm, 
                      N=100, limit=10000, seed=NA, 
                      kp=1, km=-1, hp=3, hm=-3, 
                      side="both", 
                      printSummary=TRUE) 
{
  if (!side %in% c("both", "upper", "lower")) 
    warning(gettextf("side = '%s' is not supported. Using 'both'", 
                     side), domain = NA)
  
  if(!is.na(seed))
    set.seed(seed)
  
  data <- matrix(randFunc(n=limit*N, ...), nrow=N)
  n <- ncol(data)
  
  res <- list(run=apply(data, MARGIN=1, .runLength, kp=kp, km=km, ubd=hp, lbd=hm, side=side))
  res$rls <- sapply(res$run, function(x) x$rl)
  
  res$statistic <- c(mean(res$rls), sqrt((mean(res$rls^2) - mean(res$rls))/N))
  names(res$statistic) <- c("ARL", "Std. Error")
  
  if(printSummary)
    print(res$statistic)
  
  invisible(res)
}

cusumPfaCedNorm <-  function (mean0 = 0, sd0=1, mean1=0, sd1=1, 
                              tau=10, 
                              N=100, limit=10000, seed=NA, 
                              kp=1, km=-1, hp=3, hm=-3, 
                              side="both", 
                              printSummary=TRUE) 
{
  
  
  if(!is.na(seed))
    set.seed(seed)
  
  data <- matrix(rnorm(n=(tau)*N, mean0, sd0), nrow=N)
  data <- cbind(data, matrix(rnorm(n=(limit-tau)*N, mean1, sd1), nrow=N))
  
  .cusumPfaCed(data, tau, kp, km, hp, hm, side, printSummary)
  
}

cusumPfaCedBinom <-  function (size0 = 0, prob0=1, size1=0, prob1=1, 
                               tau=10, 
                               N=100, limit=10000, seed=NA, 
                               kp=1, km=-1, hp=3, hm=-3, 
                               side="both", 
                               printSummary=TRUE) 
{
  
  
  if(!is.na(seed))
    set.seed(seed)
  
  data <- matrix(rbinom(n=(tau)*N, size0, prob0), nrow=N)
  data <- cbind(data, matrix(rnorm(n=(limit-tau)*N, size1, prob1), nrow=N))
  
  .cusumPfaCed(data, tau, kp, km, hp, hm, side, printSummary)
  
}

cusumPfaCedPois <-  function (lambda0 = 0, lambda1=1, 
                              tau=10, 
                              N=100, limit=10000, seed=NA, 
                              kp=1, km=-1, hp=3, hm=-3, 
                              side="both", 
                              printSummary=TRUE) 
{
  
  
  if(!is.na(seed))
    set.seed(seed)
  
  data <- matrix(rpois(n=(tau)*N, lambda0), nrow=N)
  data <- cbind(data, matrix(rnorm(n=(limit-tau)*N, lambda1), nrow=N))
  
  .cusumPfaCed(data, tau, kp, km, hp, hm, side, printSummary)
  
}
