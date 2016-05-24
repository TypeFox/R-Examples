doQmapSSPLIN <- function(x,fobj,...){ 
  if(!any(class(fobj)=="fitQmapSSPLIN"))
    stop("class(fobj) should be fitQmapSSPLIN")  
  UseMethod("doQmapSSPLIN")
}

doQmapSSPLIN.default <- function(x,fobj,...){
  wet <-  if(!is.null(fobj$wet.day)){
    x>=fobj$wet.day
  } else {
    rep(TRUE,length(x))
  }
  out <- rep(NA,length.out=length(x))
  out[wet] <- predict(fobj$par[[1]],x[wet])$y
  out[!wet] <- 0
  if(!is.null(fobj$wet.day))
    out[out<0] <- 0
  return(out)
}


doQmapSSPLIN.matrix <- function(x,fobj,...){
  if(ncol(x)!=length(fobj$par))
    stop("'ncol(x)' and 'nrow(fobj$par$modq)' should be eaqual\n")  
  NN <- ncol(x)
  hind <- 1:NN
  names(hind) <- colnames(x)
  hf <- list()
  class(hf) <- class(fobj)
  xx <- sapply(hind,function(i){
    hf$par <- fobj$par[i]
    hf$wet.day <- fobj$wet.day[i]
    tr <- try(doQmapSSPLIN.default(x[,i],hf,...),silent=TRUE)
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

doQmapSSPLIN.data.frame <- function(x,fobj,...){
  x <- as.matrix(x)
  x <- doQmapSSPLIN.matrix(x,fobj,...)
  x <- as.data.frame(x)
  return(x)  
}
