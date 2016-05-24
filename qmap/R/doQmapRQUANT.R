doQmapRQUANT <- function(x,fobj,...){ 
  if(!any(class(fobj)=="fitQmapRQUANT"))
    stop("class(fobj) should be fitQmapRQUANT")  
  UseMethod("doQmapRQUANT")
}

doQmapRQUANT.default <- function(x,fobj,slope.bound=c(lower=0,upper=Inf),
                               type=c("linear","linear2","tricub"),...){
  type <- match.arg(type)
  fobj$par$slope <- c(max(fobj$par$slope["lower",],slope.bound["lower"]),
                      min(fobj$par$slope["upper",],slope.bound["upper"]))
  wet <-  if(!is.null(fobj$wet.day)){
    x>=fobj$wet.day
  } else {
    rep(TRUE,length(x))
  }
  out <- rep(NA,length.out=length(x)) 
  if(type%in%c("linear","linear2")){
    out[wet] <- approx(x=fobj$par$modq[,1], y=fobj$par$fitq[,1],
                       xout=x[wet], method="linear",
                       rule=2, ties=mean)$y
    if(type=="linear2"){
      if (any(k <- x > max(fobj$par$modq)))
        out[k]  <- max(fobj$par$fitq) +
          fobj$par$slope[2]*(x[k] - max(fobj$par$modq))
      if (any(k <- x < min(fobj$par$modq)))
        out[k]  <- min(fobj$par$fitq) +
          fobj$par$slope[1]*(x[k] - min(fobj$par$modq))
    } else if(type=="linear"){
      nq <- nrow(fobj$par$modq)
      largex <- x>fobj$par$modq[nq,1]
      if(any(largex)){
        max.delta <- fobj$par$modq[nq,1] - fobj$par$fitq[nq,1]
        out[largex] <- x[largex] - max.delta
      }
    }
  } else if(type=="tricub"){
    sfun <- splinefun(x=fobj$par$modq[,1], y=fobj$par$fitq[,1],
                      method="monoH.FC")
    out[wet] <- sfun(x[wet])
  } 
  out[!wet] <- 0
  if(!is.null(fobj$wet.day))
    out[out<0] <- 0
  return(out)
}


doQmapRQUANT.matrix <- function(x,fobj,...){
  if(ncol(x)!=ncol(fobj$par$modq))
    stop("'ncol(x)' and 'nrow(fobj$par$modq)' should be eaqual\n")  
  NN <- ncol(x)
  hind <- 1:NN
  names(hind) <- colnames(x)
  hf <- list()
  class(hf) <- class(fobj)
  xx <- sapply(hind,function(i){
    ## hf <- fobj
    hf$par$modq <- matrix(fobj$par$modq[,i],ncol=1)
    hf$par$fitq <- matrix(fobj$par$fitq[,i],ncol=1)
    hf$par$slope <- matrix(fobj$par$slope[,i],ncol=1,
                           dimnames=list(c("lower","upper"),NULL))
    hf$wet.day <- fobj$wet.day[i]
    tr <- try(doQmapRQUANT.default(x[,i],hf,...),silent=TRUE)
    if(class(tr)=="try-error"){
      warning("Quantile mapping for ",names(hind)[i],
              " failed NA's produced.")
      tr <- rep(NA,nrow(x))
    }
    return(tr)
  })
  rownames(xx) <- rownames(x)
  return(xx)
}

doQmapRQUANT.data.frame <- function(x,fobj,...){
  x <- as.matrix(x)
  x <- doQmapRQUANT.matrix(x,fobj,...)
  x <- as.data.frame(x)
  return(x)  
}
