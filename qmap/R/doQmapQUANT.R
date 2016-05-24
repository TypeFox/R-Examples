doQmapQUANT <- function(x,fobj,...){ 
  if(!any(class(fobj)=="fitQmapQUANT"))
    stop("class(fobj) should be fitQmapQUANT")  
  UseMethod("doQmapQUANT")
}

doQmapQUANT.default <- function(x,fobj,type=c("linear","tricub"),...){
  type <- match.arg(type)
  wet <-  if(!is.null(fobj$wet.day)){
    x>=fobj$wet.day
  } else {
    rep(TRUE,length(x))
  }
  out <- rep(NA,length.out=length(x))
  if(type=="linear"){
    out[wet] <- approx(x=fobj$par$modq[,1], y=fobj$par$fitq[,1],
                       xout=x[wet], method="linear",
                       rule=2, ties=mean)$y
    nq <- nrow(fobj$par$modq)
    largex <- x>fobj$par$modq[nq,1]
    if(any(largex)){         
      max.delta <- fobj$par$modq[nq,1] - fobj$par$fitq[nq,1]
      out[largex] <- x[largex] - max.delta
    }
  } else if(type=="tricub"){
    sfun <- splinefun(x=fobj$par$modq[,1], y=fobj$par$fitq[,1],
                      method="monoH.FC")#,method="monoH.FC")   
    out[wet] <- sfun(x[wet])
  }
  out[!wet] <- 0
  if(!is.null(fobj$wet.day))
    out[out<0] <- 0
  return(out)
}


doQmapQUANT.matrix <- function(x,fobj,...){
  if(ncol(x)!=ncol(fobj$par$modq))
    stop("'ncol(x)' and 'nrow(fobj$par$modq)' should be eaqual\n")  
  NN <- ncol(x)
  hind <- 1:NN
  names(hind) <- colnames(x)
  hf <- list()
  class(hf) <- class(fobj)
  xx <- sapply(hind,function(i){
    hf$par$modq <- matrix(fobj$par$modq[,i],ncol=1)
    hf$par$fitq <- matrix(fobj$par$fitq[,i],ncol=1)
    hf$wet.day <- fobj$wet.day[i]
    tr <- try(doQmapQUANT.default(x[,i],hf,...),silent=TRUE)
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

doQmapQUANT.data.frame <- function(x,fobj,...){
  x <- as.matrix(x)
  x <- doQmapQUANT.matrix(x,fobj,...)
  x <- as.data.frame(x)
  return(x)  
}
