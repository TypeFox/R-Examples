################################################################
## Computing properties of CUSUM charts with known parameters ##
################################################################

hitprob_sim <- function(c,nrep=1000,n=1e3,robs){
  mean(replicate(nrep,{
    R <- cumsum(robs(n))
    S <- R-cummin(R)
    any(S>=c)
  }))
}

calibratehitprob_sim <- function(hprob=0.01,nrep=1000,n=1e3,robs){
  x <- (replicate(nrep,{
    R <- cumsum(robs(n))
    S <- R-cummin(R)
    max(S)
  }))
  quantile(x,1-hprob)
}



ARL_sim <- function(c,nrep=1000,maxsteps=1e4,robs){
  mean(replicate(nrep,{
    R <- cumsum(robs(maxsteps))
    S <- R-cummin(R)
    w <- S>=c
    if (any(w))
      return(min(which(w)))
    else
      return(maxsteps*1.5)
  }))
}

getQ <- function(c,gridpoints,pobs){
  p <- c(0,(0.5+(0:(gridpoints-2)))*c/(gridpoints-1))
  ptarget <- p+c/(gridpoints-1)/2
  ptarget[1] <- ptarget[1]- c/(gridpoints-1)/4
  sapply(p,function(x) {res <- pobs(ptarget-x); c(res[1],diff(res))})
}

ARL_CUSUM_Markovapprox <- function(c,gridpoints=100,pobs,gridpointsmax=500){
    if (c<=0) return(1);
    if (pobs(0)>=1) return(Inf);
    gridpointsakt <- gridpoints
    while(gridpointsakt<gridpoints*8 && gridpointsakt<gridpointsmax){
        gridsize <- c/(gridpointsakt-1)/2 
        if (pobs(gridsize)<0.99) break;
        gridpointsakt <- min(gridpointsmax,gridpointsakt*2)
        if (gridpointsakt==gridpointsmax) {
            warning(paste("Adjustment of grid size reached maximum of",gridpointsmax, "grid points."))
        }
    }
    if (pobs(c/(gridpointsakt-1)/2.01)>=1) return(Inf); ##cannot reach third grid point
    Q <- getQ(c,gridpointsakt,pobs)
    tryCatch(rep(1,gridpointsakt)%*%solve(diag(rep(1,gridpointsakt))-Q,c(1,rep(0,gridpointsakt-1))),error=function(e) Inf)
}


calibrateARL_Markovapprox<- function(ARL=1000,f=ARL_CUSUM_Markovapprox,pobs,...){
  cmax <- 4;
  while(f(pobs=pobs,cmax,...)<ARL) cmax <- cmax*1.5
  cmin <- cmax/1.5;
  if (cmin<4){
      while(f(pobs=pobs,cmin,...)>ARL&&cmin>1e-10) cmin <- cmin/2
  }
  if (cmin<=1e-10){
    stop("Calibration to an ARL of ",ARL," not possible.",
         "Hitting probability with threshold ", format(cmin), " gives an ARL of ",
         format(f(pobs=pobs,cmin,...)),".")
  }
  warningslist <- c();
  withCallingHandlers({
      resuniroot <- uniroot(function(x) {      
          res <- f(pobs=pobs,c=x,...)-ARL;
          ifelse(res==Inf, .Machine$double.xmax, res) ## to avoid warning messages replacing Inf by largest double
      },lower=cmin,upper=cmax)
  },warning=function(w){
      warningslist <<- c(warningslist,w$message)
      invokeRestart("muffleWarning");
  })
  if (length(warningslist)>0){
      tab <- table(warningslist)
      for (i in 1:length(tab)){
         if (tab[[i]]>1)
             warning(paste(tab[[i]],"times:",names(tab)[i]))
         else
             warning(names(tab)[i])
     }
  }
  resuniroot$root
}


matrix.power <- function(mat, n)
{
  if (n == 1) return(mat)
  result <- diag(1, ncol(mat))
  while (n > 0) {
    if (n %% 2 != 0) {
      result <- result %*% mat
      n <- n - 1
    }
    mat <- mat %*% mat
    n <- n / 2
  }
  return(result)
}


hitprob_CUSUM_Markovapprox <- function(pobs,c,n,gridpoints=100,gridpointsmax=200){
    if (c<=0) return(1);
    if (pobs(0)>=1) return(0);
    gridpointsakt <- gridpoints
    while(gridpointsakt<gridpoints*8 && gridpointsakt<gridpointsmax){
        gridsize <- c/(gridpointsakt-1)/2 
        if (pobs(gridsize)<0.99) break;
        gridpointsakt <- min(gridpointsmax,gridpointsakt*2)
        if (gridpointsakt==gridpointsmax) warning(paste("Adjustment of grid size reached maximum of",gridpointsmax, "grid points."))
    }
    if (pobs(c/(gridpointsakt-1)/2.01)>=1) return(Inf); ##cannot reach third grid point
  Q <- getQ(c,gridpoints=gridpointsakt,pobs=pobs)
  Q <- matrix.power(Q,n)
  1-rep(1,gridpointsakt)%*%Q%*%c(1,rep(0,gridpointsakt-1))
}


calibratehitprob_Markovapprox<- function(hprob=0.01,n,f=hitprob_CUSUM_Markovapprox,pobs,...){
  cmax <- 1;
  while(f(pobs=pobs,n=n,cmax,...)>hprob) cmax <- cmax*2
  cmin <- cmax/2;
  if (cmin<1){
    while(f(pobs=pobs,n=n,cmin,...)<hprob&&cmin>1e-10) cmin <- cmin/2
  }
  if (cmin<=1e-10){
    stop("Calibration to a hitting probability of ",hprob,
         " not possible. Hitting probability with threshold ",format(cmin),
         " is ",format(f(pobs=pobs,n=n,cmin,...)),".")
  }
  uniroot(function(x) f(pobs=pobs,n=n,c=x,...)-hprob,lower=cmin,upper=cmax)$root
}


