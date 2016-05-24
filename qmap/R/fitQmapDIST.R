fitQmapDIST <- function(obs,mod,...)
  UseMethod("fitQmapDIST")

fitQmapDIST.default <- function(obs,mod,distr="berngamma",start.fun,
                                qstep=NULL,mlepar,...){ 
### identification of the parameters for parametric
### quantile mapping using the function 'mledist' from
### package 'fitdistrplus'.
###
### obs: oserved time series
### mod: modelled time series
### distr: A character string '"name"' naming a distribution for which
###        the corresponding density function 'dname' and the
###        corresponding distribution 'pname' must be classically
###          defined. 
### start.fun: function estimating starting values for
###            parameter optimisation. Default starting
###            values are provided for "mixgamma" and the
###            distributions mentioned in the documentation
###            of 'mledist'.
### qstep: mod and obs are aggregated to quantiles
###        before model identification as:
###        quantile(x,probs=seq(0,1,by=qstep).
### mlepar: further parameters passed to 'mledist' note that
###         'start' may be overwritten by 'start.fun'
### ... : further arguments passed to 'mledist' and 'optim'.
  obs <- na.omit(obs)
  mod <- na.omit(mod)
  if(!is.null(qstep)){
    if(qstep<1&qstep>0){
      obs <- quantile(obs,probs=seq(0,1,by=qstep),type=8)
      mod <- quantile(mod,probs=seq(0,1,by=qstep),type=8)
    } else{
      stop("'qstep' shoub be NULL or in the 'qstep < 1 & qstep > 0' intervall")
    }
  }
  start.fun <- if(missing(start.fun)){
    sf <- paste("start",distr,sep="")
    if(exists(sf,mode="function"))
      match.fun(sf) else function(x)return(NULL)
  } else start.fun
  if(missing(mlepar))
    mlepar <- list() 
  if(!any(names(mlepar)%in%c("lower","upper"))){
    if(distr=="berngamma"){
      mm <- .Machine$double.eps
      mlepar$lower=c(0,mm,mm)
      mlepar$upper=c(1,Inf,Inf)
    }
  }
  mlepar$distr <- distr
  mlepar$start <- start.fun(obs)
  ow <- options(warn=-1)
  obs.fit <- do.call(mledist,c(list(obs),mlepar))
  options(ow)
    if(obs.fit$convergence > 0) 
      stop("'mledist' failed to estimate parameters for 'obs' with the error code ",
           obs.fit$convergence, "\n")
  mlepar$start <- start.fun(mod)
  ow <- options(warn=-1)
  mod.fit <- do.call(mledist,c(list(mod),mlepar))
  options(ow)
  if(mod.fit$convergence > 0) 
      stop("'mledist' failed to estimate parameters for 'mod' with the error code ",
           mod.fit$convergence, "\n")
  obs.fit <- obs.fit$estimate
  mod.fit <- mod.fit$estimate
  pdistnam <- paste("p",distr,sep="")
  pdist <- match.fun(pdistnam)
  parg <- names(formals(pdist))
  qdistnam <- paste("q",distr,sep="")
  qdist <- match.fun(qdistnam)
  qarg <- names(formals(qdist))  
  nm.o <- names(obs.fit)
  nm.m <- names(mod.fit)
  m.o <- match(nm.o,qarg)
  m.m <- match(nm.m,parg)
  arg.o <- paste(qarg[m.o],"o",sep=".")
  arg.m <- paste(parg[m.m],"m",sep=".")
  ## build transfer funciton
  ## arguments...
  tfun.arg <- paste(c("x",arg.o,arg.m),collapse=",")
  ## body...
  fl1 <- paste(paste(parg[m.m],arg.m,sep="="),collapse=",")
  fl1 <- paste("pp<-",pdistnam,"(x,",fl1,")")
  fl1 <- parse(text=fl1)
  fl2 <- paste(paste(qarg[m.o],arg.o,sep="="),collapse=",")
  fl2 <- paste("qq<-",qdistnam,"(pp,",fl2,")")
  fl2 <- parse(text=fl2)
  tfun.body <- as.call(c(as.name("{"),
                         fl1,
                         fl2,
                         expression(return(qq))))
  ## put the funciton together
  tfun <- paste("function(",tfun.arg,"){}")
  tfun <- eval(parse(text=tfun))
  body(tfun) <- tfun.body
  ## prepare output
  ppar <- c(obs.fit,mod.fit)
  ppar <- matrix(ppar,nrow=1)
  colnames(ppar) <- c(arg.o,arg.m)
  op <- list(tfun=tfun,
             par=ppar)
  class(op) <- "fitQmapDIST"
  return(op)
}


fitQmapDIST.matrix <- function(obs,mod,...){
  if(ncol(mod)!=ncol(obs))
    stop("'mod' and 'obs' need the same number of columns")
  NN <- ncol(mod)
  hind <- 1:NN
  names(hind) <- colnames(mod)
  xx <- lapply(hind,function(i){
    tr <- try(fitQmapDIST.default(obs=obs[,i],mod=mod[,i],...),
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
  xx <- list(tfun=tfun,
             par=ppar)
  class(xx) <- "fitQmapDIST"
  return(xx)    
}

fitQmapDIST.data.frame <- function(obs,mod,...){
  obs <- as.matrix(obs)
  mod <- as.matrix(mod)
  fitQmapDIST.matrix(obs,mod,...)
}
