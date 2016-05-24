fitQmapPTF <- function(obs,mod,...)
  UseMethod("fitQmapPTF")


fitQmapPTF <- function(obs,mod,...)
  UseMethod("fitQmapPTF")

fitQmapPTF.default <- function(obs,mod,
                               transfun=c("power","linear","expasympt",
                                 "scale","power.x0","expasympt.x0"),
                               wet.day=TRUE,cost=c("RSS","MAE"),
                               qstep=0.001,opar,...){
### identify parameters of transfere functions for
### bias-correction of RCM
###
### obs, mod: vector with observed/modelled values
### transfun:
### -- A -- character string indicating transfere functuion
### --> power: y=b*x^c
### --> linear: y = a + b*x
### --> scale: y = b*x
### --> expasympt: y=(a+b*x)*(1-exp(-x/tau))
### --> mixgamma: "quantile mapping" with mixed
###     gamma distributuib
### -- B -- an function with 'x' as first argumment
### --> e.g.: funciton(x,a,b){...}
### wet.day: threshold for minimal precipitation for 'wett day'
### cost: cost function to minimize.
### --> RSS: residual sum of squares
### --> MAE: median absolute error
### qstep: mod and obs are aggregated to quantiles
###        before model identification as:
###        quantile(x,probs=seq(0,1,by=qstep).
### opar: a list with arguments passed to 'optim'. Note that 'method'
###      is chosen automatically. The 'par' argument (starting values)
###      has defaults for the predefined funcitons
  cost <- match.arg(cost)
  if(missing(opar))
    opar <- list()
  obs <- na.omit(obs)
  mod <- na.omit(mod)
  if(is.null(qstep)){
    if(length(obs)!=length(mod)){
      hn <- min(length(obs),length(mod))
      obsq <- quantile(obs,seq(0,1,length.out=hn),type=8)
      modq <- quantile(mod,seq(0,1,length.out=hn),type=8)
    } else {
      obsq <- sort(obs)
      modq <- sort(mod)
    }
  } else if(qstep<1&qstep>0) {
    obsq <- quantile(obs,probs=seq(0,1,by=qstep),type=8)
    modq <- quantile(mod,probs=seq(0,1,by=qstep),type=8)
  } else {
    stop("'qstep' shoub be NULL or in the 'qstep < 1 & qstep > 0' intervall")
  }
  if(is.numeric(wet.day)){
    q0 <- obsq>=wet.day
    obsq <- obsq[q0]
    modq <- modq[q0]
  } else if(is.logical(wet.day)){
    if(wet.day){
      q0 <- obsq>0
      obsq <- obsq[q0]
      modq <- modq[q0]
      wet.day <- modq[1]
      names(wet.day) <- NULL
    } else {
      wet.day <- NULL
    }
  } else {
    stop("'wet.day' should be 'numeric' or 'logical'")
  } 
  tfun <- if(is.character(transfun)){
    transfun <- match.arg(transfun)   
    switch(transfun,
           power.x0=function(x,a,b,x0){
             a*(x-x0)^b
           },
           power=function(x,a,b){
             a*x^b
           },
           linear=function(x,a,b){
             x <- a+b*x
           },
           expasympt.x0=function(x,a,b,x0,tau){
             x <- (a+b*x)*(1-exp(-(x-x0)/tau))
             return(x)
           },
           expasympt=function(x,a,b,tau){
             gz <- x>0
             x[!gz] <- 0
             x[gz] <- (a+b*x[gz])*(1-exp(-(x[gz])/tau))
             return(x)
           },
           scale=function(x,b){
             x*b
           })
  } else if(is.function(transfun)){
    if(is.null(opar$par))
      stop("'opar$par' needs to be supplied if 'transfun' is a function\n")
    fpp <- formals(transfun)
    fpp <- names(fpp)
    if(fpp[1]!="x")
      stop("first argument of 'transfun' should be 'x'\n")
    if(!all(fpp[-1]%in%names(opar$par)))
      stop("names of 'opar$par' should correspond to arguments of 'transfun'\n")
    transfun
  } else {
    stop("'transfun' needs to be character or function")
  }
  if(is.null(opar$par)){
    opar$par <- switch(transfun,
                       power.x0=c(a=1,b=1,x0=0),
                       power=c(a=1,b=1),
                       linear={
                         pp <- coef(lm(obsq~modq))
                         names(pp) <- c("a","b")
                         pp},
                       expasympt.x0={
                         pp <- coef(lm(obsq~modq))
                         pp <- c(pp,0,0)
                         names(pp) <- c("a","b","x0","tau")
                         pp
                       },
                       expasympt={
                         pp <- coef(lm(obsq~modq))
                         pp <- c(pp,0)
                         names(pp) <- c("a","b","tau")
                         pp
                       },
                       scale=c(b=mean(obsq)/mean(modq))
                       )
  }
  mfun <- switch(cost,
                 RSS=function(par,x,y){
                   rr <- y - do.call(tfun,c(list(x),as.list(par)))
                   sum(rr^2)
                 },
                 MAE=function(par,x,y){
                   rr <- y - do.call(tfun,c(list(x),as.list(par)))
                   mean(abs(rr))
                 })
  opar$fn <- mfun
  opar$x <- modq
  opar$y <- obsq
  if(any(names(opar)%in%c("lower","upper"))){
    opar$method <- "L-BFGS-B"
  } else{
    if(length(opar$par)==1){
      opar$method <- "BFGS"
    } else {
      opar$method <- "Nelder-Mead"
    }
  }
  opt <- try(do.call(optim,opar),silent=TRUE)
  if(class(opt)=="try-error"){
    warning("optim 'method'",opar$method,
            "failed. Optimizing with 'SANN'\n",
            "possibly unstable result")
    opar$method <- "SANN"
    opar$lower <- NULL
    opar$upper <- NULL
    opt <- do.call(optim,opar)
  } 
  ppar <- t(as.matrix(opt$par))
  op <- list(tfun=tfun,
             par=ppar,
             wet.day=wet.day
             )
  class(op) <- "fitQmapPTF"
  return(op)  
}

fitQmapPTF.matrix <- function(obs,mod,...){
  if(ncol(mod)!=ncol(obs))
    stop("'mod' and 'obs' need the same number of columns")
  NN <- ncol(mod)
  hind <- 1:NN
  names(hind) <- colnames(mod)
  xx <- lapply(hind,function(i){
    tr <- try(fitQmapPTF.default(obs=obs[,i],mod=mod[,i],...),
              silent=TRUE)
    if(any(class(tr)=="try-error")){
      warning("model identification for ",names(hind)[i],
              " failed\n NA's produced")
      NULL
    } else{
      tr
    }
  })
  xx.NULL <- sapply(xx,is.null)
  tfun <- xx[!xx.NULL][[1]]$tfun
  ppar <- lapply(xx,function(x)x$par)
  parnam <- colnames(ppar[!xx.NULL][[1]])
  npar <- ncol(xx[!xx.NULL][[1]]$par)
  npar <- matrix(NA,ncol=npar,
                 dimnames=list(NULL,parnam))
  ppar[xx.NULL] <- list(npar)
  ppar <- do.call(rbind,ppar)
  rownames(ppar) <- names(xx)
  wday <- lapply(xx,function(x)x$wet.day)
  wday[xx.NULL] <- NA
  wday <- do.call(c,wday)
  xx <- list(tfun=tfun,
             par=ppar,
             wet.day=wday)
  class(xx) <- "fitQmapPTF"
  return(xx)    
}

fitQmapPTF.data.frame <- function(obs,mod,...){
  obs <- as.matrix(obs)
  mod <- as.matrix(mod)
  fitQmapPTF.matrix(obs,mod,...)
}


fitQmapPTF.matrix <- function(obs,mod,...){
  if(ncol(mod)!=ncol(obs))
    stop("'mod' and 'obs' need the same number of columns")
  NN <- ncol(mod)
  hind <- 1:NN
  names(hind) <- colnames(mod)
  xx <- lapply(hind,function(i){
    tr <- try(fitQmapPTF.default(obs=obs[,i],mod=mod[,i],...),
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
  tfun <- xx[!xx.NULL][[1]]$tfun
  ppar <- lapply(xx,function(x)x$par)
  parnam <- colnames(ppar[!xx.NULL][[1]])
  npar <- ncol(xx[!xx.NULL][[1]]$par)
  npar <- matrix(NA,ncol=npar,
                 dimnames=list(NULL,parnam))
  ppar[xx.NULL] <- list(npar)
  ppar <- do.call(rbind,ppar)
  rownames(ppar) <- names(xx)
  wday <- lapply(xx,function(x)x$wet.day)
  wday[xx.NULL] <- NA
  wday <- do.call(c,wday)
  xx <- list(tfun=tfun,
             par=ppar,
             wet.day=wday)
  class(xx) <- "fitQmapPTF"
  return(xx)    
}

fitQmapPTF.data.frame <- function(obs,mod,...){
  obs <- as.matrix(obs)
  mod <- as.matrix(mod)
  fitQmapPTF.matrix(obs,mod,...)
}
