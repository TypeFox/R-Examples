dpoilog <- function(n, mu, sig, trunc=FALSE) {
  if (length(mu)>1 | length(sig)>1) stop('vectorization of mu and sig is not supported') 
  if (any((n[n!=0]/trunc(n[n!=0]))!=1)) stop('all n must be integers')
  if (!all(is.finite(c(mu,sig)))) stop('all parameters should be finite')
  if (sig<=0) stop('sig is not larger than 0')
  if (trunc & all(n == 0)) warning('all n will be truncated')
  
  pdfs <- .dpoilog(n, mu, sig)
  if (trunc) {
    pdfs <- pdfs / (1.0 - .dpoilog(0, mu, sig))
    pdfs[n == 0] <- 0.0
  }
  pdfs
}

.dpoilog <- function(n, mu, sig) {
  .C('poilog', as.integer(n), as.double(mu), as.double(sig),
     n_obs=as.integer(length(n)), val=double(length(n)))$val
}

rpoilog <- function(S, mu, sig, condS=FALSE, keep0=FALSE){
   sim <- function(nr){
     lamx <- rnorm(nr)
     x <- rpois(nr,exp(sig*lamx+mu))
     if (!keep0) x <- x[x>0]
     return(x)
   }
   
   if (S<1) stop('S is not positive')
   if (!is.finite(S)) stop('S is not finite')
   if ((S/trunc(S))!=1) stop('S is not an integer')
   if (sig<0) stop('sig is not positive')
   
   if (condS) {
     simVec <- vector('numeric',0)
     fac <- 2
     nr  <- S
     while (length(simVec)<S){
       simvals <- sim(nr*fac)
       simVec <- c(simVec,simvals)
       fac <- (1/(length(simvals)/(nr*fac)))*2
       fac <- ifelse(is.finite(fac),fac,1000)
       nr <- S-length(simvals)
     }
     simVec <- simVec[1:S]
   }
   
   else simVec <- sim(S)
   return(simVec)
}

poilogMLE <- function(n, nboot=0, trunc=TRUE, method='L-BFGS-B', start.mu=-1.0, start.sig=1.0,
  control=list(fnscale=length(n)), ...) {

  if (is.matrix(n) | (is.data.frame(n))) {
    stop(paste('n has',ncol(n),'colums, supply a vector or use function bipoilogMLE',sep=' ')) 
  }
  
  # truncate the input
  if (trunc) n <- n[n > 0]
  
  # guess innocuous start values
  startVals=c(mu=start.mu, sig=start.sig)

  # dereplicate
  un <- unique(n)
  nr <- rep(NA,length(un))
  for (i in 1:length(un)){ nr[i] <- sum(n%in%un[i]) }
  
  # negative log likelihood is the objective
  lnL <- function(z) {
    -sum((log(dpoilog(un, z[1], exp(z[2]), trunc=trunc)))*nr)
  }
  fit <- optim(startVals, lnL, control=control, method=method, lower=-20, upper=20, ...)
  
  if (fit$convergence!=0){
    if (fit$convergence==1) stop('the iteration limit has been reached!   try different startVals or increase maxit') 
    if (fit$convergence==10) stop('degeneracy of the Nelder Mead simplex ....')
    else stop(paste('unknown error in optimization', fit$message))
  } 
  
  fit$par <- c(as.numeric(fit$par),1-dpoilog(0,fit$par[1],exp(fit$par[2])))
  est <- list('par'=c('mu'=fit$par[1],'sig'=exp(fit$par[2])),'p'=fit$par[3],'logLval'=-fit$value,'gof'=NULL,boot=NULL)
  
  if (nboot>0){
    cat(paste('estimates: mu: ',round(fit$par[1],3),' sig ',round(exp(fit$par[2]),3),sep=''),'\n')
    cat('********     bootstrapping    ********\n')
    bMat <- matrix(NA,nboot,3)
    colnames(bMat) <- c('mu','sig2','logLval')
    count <- 0
    kat  <- seq(0,nboot,by=100)
    bStartVals <- fit$par
    while (count<nboot){
      bfit <- un <- nr <- NA
      # simulations are conditonal on the number of species in the observed data set :
      sim <- rpoilog(length(n),fit$par[1],exp(fit$par[2]),condS=TRUE,keep0=!trunc)
      un <- unique(sim)
      nr <- rep(NA,length(un))
      for (i in 1:length(un)){ nr[i] <- sum(sim%in%un[i]) }
      bfit <- try(optim(bStartVals,lnL,nr=nr,control=control,method=method),silent=TRUE)
      if (class(bfit)!='try-error'){
        count <- count+1
        bMat[count,] <- c(bfit$par[1],exp(bfit$par[2]),-bfit$value)
        if (count%in%kat) cat('   boot',count,'of',nboot,'\n')
      }
    }
    est$boot <- data.frame(bMat)
    est$gof  <- which(sort(c(est$logLval,bMat[,3]))==est$logLval)/nboot
  }  
  return(est)   
}
