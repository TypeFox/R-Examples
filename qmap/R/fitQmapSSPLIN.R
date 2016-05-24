fitQmapSSPLIN <- function(obs,mod,...)
  UseMethod("fitQmapSSPLIN")

fitQmapSSPLIN.default <- function(obs,mod,wet.day=TRUE,qstep=0.01,spline.par,...){
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
  xs <- if(!is.null(qstep)){
    quantile(xs, probs=seq(0,1,by=qstep),type=8) 
  } else sort(xs)
  ys <- if(!is.null(qstep)){
    quantile(ys, probs=seq(0,1,by=qstep),type=8)    
  } else sort(ys)
  if(missing(spline.par))
    spline.par <- list()
  spline.par$x <- xs
  spline.par$y <- ys
  fit <- do.call(smooth.spline,spline.par)
  ## 13.03.2012
  ## force spline to be monotonic:
  fit$fit$coef <- cummax(fit$fit$coef)
  op <- list(par=list(fit$fit),
             wet.day=wet.day)
  class(op) <- "fitQmapSSPLIN"
  return(op)
}

fitQmapSSPLIN.matrix <- function(obs,mod,...){
  if(ncol(mod)!=ncol(obs))
    stop("'mod' and 'obs' need the same number of columns") 
  NN <- ncol(mod)
  hind <- 1:NN
  names(hind) <- colnames(mod)
  xx <- lapply(hind,function(i){
    tr <- try(fitQmapSSPLIN.default(obs=obs[,i],mod=mod[,i],...),
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
  ppar <- lapply(xx,function(x)x$par[[1]])
  wday <- lapply(xx,function(x)x$wet.day)
  wday[xx.NULL] <- NA
  wday <- do.call(c,wday)
  xx <- list(par=ppar,
             wet.day=wday)
  class(xx) <- c("fitQmapSSPLIN")
  return(xx)    
}


fitQmapSSPLIN.data.frame <- function(obs,mod,...){
  obs <- as.matrix(obs)
  mod <- as.matrix(mod)
  fitQmapSSPLIN.matrix(obs,mod,...)
}
