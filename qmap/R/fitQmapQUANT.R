fitQmapQUANT <- function(obs,mod,...)
  UseMethod("fitQmapQUANT")

fitQmapQUANT.default <- function(obs,mod,wet.day=TRUE,qstep=0.01,
                                 nboot = 1,...){
### fits an nonparametric transfere function for quantile
### mapping using local linear least square regression
###
### Based on code from John Bjornar Bremnes
  ys <- na.omit(obs)
  xs <- na.omit(mod)
  if(length(xs)!=length(ys)){
    hn <- min(length(xs),length(ys))
    ## quantile algorithm 'type=8' appeares to
    ## be reccomended. See help(quantile) and
    ## Hyndman & Fan (1996) mentioned therein
    ys <- quantile(ys,seq(0,1,length.out=hn),type=8)
    xs <- quantile(xs,seq(0,1,length.out=hn),type=8)
  } else {
    xs <- sort(xs)
    ys <- sort(ys)
  } 
  if(is.numeric(wet.day)){
    q0 <- ys>=wet.day
    ys <- ys[q0]
    xs <- xs[q0]
  } else if(is.logical(wet.day)){
    if(wet.day){
      q0 <- ys>0
      ys <- ys[q0]
      xs <- xs[q0]
      wet.day <- xs[1]
      names(wet.day) <- NULL
    } else {
      wet.day <- NULL
    }
  } else {
    stop("'wet.day' should be 'numeric' or 'logical'")
  } 
  nn <- length(ys)  
  ## define predictor values at which a fit is wanted
  newx <- quantile(xs, probs=seq(0,1,by=qstep),type=8)
  fit   <- array(NA, dim=c(length(newx),2,nboot))
  if(nboot > 1){
    booty <- replicate(nboot,sample(ys,size=nn,replace=FALSE))
    booty <- apply(booty,2,quantile,
                   probs=seq(0,1,by=qstep),type=8)
    booty <- rowMeans(booty,na.rm=TRUE)
  } else {
    booty <- quantile(ys,probs=seq(0,1,by=qstep),type=8)
  }
  ppar <- list(modq=matrix(newx,ncol=1),
               fitq=matrix(booty,ncol=1)) 
  op <- list(par=ppar,
             wet.day=wet.day)
  class(op) <- "fitQmapQUANT"
  return(op)
}

fitQmapQUANT.matrix <- function(obs,mod,...){
  if(ncol(mod)!=ncol(obs))
    stop("'mod' and 'obs' need the same number of columns") 
  NN <- ncol(mod)
  hind <- 1:NN
  names(hind) <- colnames(mod)
  xx <- lapply(hind,function(i){
    tr <- try(fitQmapQUANT.default(obs=obs[,i],mod=mod[,i],...),
              silent=TRUE)
    if(any(class(tr)=="try-error")){
      warning("model identification for ",names(hind)[i],
              " failed\n NA's produced.")
      NULL
    } else{
      tr
    }
  })
  xx.NULL <- sapply(xx,is.null) 
  modq <- lapply(xx,function(x)x$par$modq)
  fitq <- lapply(xx,function(x)x$par$fitq)
  nq <- nrow(xx[!xx.NULL][[1]]$par$modq)
  nq <- matrix(NA,nrow=nq)
  modq[xx.NULL] <- list(nq)
  fitq[xx.NULL] <- list(nq)
  modq <- do.call(cbind,modq)
  fitq <- do.call(cbind,fitq)
  colnames(modq) <- names(xx)
  colnames(fitq) <- names(xx)
  wday <- lapply(xx,function(x)x$wet.day)
  wday[xx.NULL] <- NA
  wday <- do.call(c,wday)
  xx <- list(par=list(modq=modq,fitq=fitq),
             wet.day=wday)
  class(xx) <- c("fitQmapQUANT")
  return(xx)    
}


fitQmapQUANT.data.frame <- function(obs,mod,...){
  obs <- as.matrix(obs)
  mod <- as.matrix(mod)
  fitQmapQUANT.matrix(obs,mod,...)
}
