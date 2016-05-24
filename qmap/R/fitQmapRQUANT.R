fitQmapRQUANT <- function(obs,mod,...)
  UseMethod("fitQmapRQUANT")

fitQmapRQUANT.default <- function(obs,mod,wet.day=TRUE,qstep=0.01,
                                nlls = 10,nboot = 10,...){
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
  nlls2 <- min(nlls, nn)
  for (j in 1:nboot) {
    if(nboot==1){
      xss <- sort.default(xs)
      yss <- sort.default(ys)
    } else {
      xss <- sort.default(sample(xs, size=nn, replace=FALSE))
      yss <- sort.default(sample(ys, size=nn, replace=FALSE))
    }
    for (i in 1:length(newx)) {
      xc     <- xss - newx[i]
      mdist  <- sort.default(abs(xc))[nlls2]
      k      <- abs(xc) <= mdist
      xc     <- cbind(1, xc[k])
      a      <- t.default(xc) %*% xc
      if ( abs(d <- a[1,1]*a[2,2] - a[1,2]*a[2,1]) < 1e-10 )  # singular? <=> det(x^tx)=0?
        fit[i,1,j] <- mean(yss[k])
      else {
        b         <- t.default(xc) %*% yss[k]
        fit[i,,j] <- c(b[1]*a[2,2]-b[2]*a[1,2], b[2]*a[1,1]-b[1]*a[2,1]) / d
      }
    }
  }
  if (!is.null(wet.day))
    if (any(k <- fit[,1,] < 0))
      fit[,1,][which(k)] <- 0
  fitted <- if(nboot>1){
    rowMeans(fit[,1,],na.rm=TRUE)
  } else fit[,1,] 
  slope <- if(nboot>1){
    rowMeans(fit[,2,][,,drop=FALSE],na.rm=TRUE)
  } else fit[,2,]
  slope <- slope[range(which(!is.na(slope)))]
  slope <- matrix(slope,ncol=1,
                  dimnames=list(c("lower","upper"),NULL))
  ppar <- list(modq=matrix(newx,ncol=1),
               fitq=matrix(fitted,ncol=1),
               slope=slope) 
  op <- list(par=ppar,
             wet.day=wet.day)
  class(op) <- "fitQmapRQUANT"
  return(op)
}

fitQmapRQUANT.matrix <- function(obs,mod,...){
  if(ncol(mod)!=ncol(obs))
    stop("'mod' and 'obs' need the same number of columns") 
  NN <- ncol(mod)
  hind <- 1:NN
  names(hind) <- colnames(mod)
  xx <- lapply(hind,function(i){
    tr <- try(fitQmapRQUANT.default(obs=obs[,i],mod=mod[,i],...),
              silent=TRUE)
    if(any(class(tr)=="try-error")){
      warning("model identification for ",names(hind)[i],
              " failed. NA's produced")
      NULL
    } else{
      tr
    }
  })
  xx.NULL <- sapply(xx,is.null) 
  ## tfun <- xx[!xx.NULL][[1]]$tfun
  modq <- lapply(xx,function(x)x$par$modq)
  fitq <- lapply(xx,function(x)x$par$fitq)
  slope <- lapply(xx,function(x)x$par$slope)
  nq <- nrow(xx[!xx.NULL][[1]]$par$modq)
  nq <- matrix(NA,nrow=nq)
  modq[xx.NULL] <- list(nq)
  fitq[xx.NULL] <- list(nq)
  slope[xx.NULL] <- list(matrix(NA,ncol=1,nrow=2,
                           dimnames=list(c("lower","upper"),NULL)))
  modq <- do.call(cbind,modq)
  fitq <- do.call(cbind,fitq)
  slope <- do.call(cbind,slope)
  colnames(modq) <- names(xx)
  colnames(fitq) <- names(xx)
  colnames(slope) <- names(xx)
  wday <- lapply(xx,function(x)x$wet.day)
  wday[xx.NULL] <- NA
  wday <- do.call(c,wday)
  xx <- list(par=list(modq=modq,fitq=fitq,slope=slope),
             wet.day=wday)
  class(xx) <- c("fitQmapRQUANT")
  return(xx)    
}


fitQmapRQUANT.data.frame <- function(obs,mod,...){
  obs <- as.matrix(obs)
  mod <- as.matrix(mod)
  fitQmapRQUANT.matrix(obs,mod,...)
}
