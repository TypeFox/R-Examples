iim <-
function(x, which.ones=NULL)
{
  if(NROW(x)==1){
    n <- x
    mIIS <- matrix(0,n,n)
    diag(mIIS) <- 1
    colnames(mIIS) <- paste("iis", 1:n, sep="")
    if(!is.null(which.ones)){ mIIS <- mIIS[,which.ones] }
    mIIS <- as.zoo(mIIS)
  }else{
    n <- NROW(x)
    mIIS <- matrix(0,n,n)
    diag(mIIS) <- 1
    x <- as.zoo(x)
    x.index <- index(x)
    xIsRegular <- is.regular(x, strict=TRUE)
    if(xIsRegular){
      xIndexObs <- floor(as.numeric(x.index))
      xCycle <- as.numeric(cycle(x))
      xIndexAsChar <- paste(xIndexObs, "(", xCycle, ")", sep="")
      xFrequency <- frequency(x)
    }else{
      xIndexAsChar <- as.character(x.index)
    }
    colnames(mIIS) <- paste("iis",
      xIndexAsChar, sep="")
    mIIS <- zoo(mIIS, order.by=x.index)
    if(xIsRegular){ mIIS <- as.zooreg(mIIS) }
    if(!is.null(which.ones)){
      where.indicators <- which(index(mIIS) %in% which.ones)
      if(length(where.indicators > 0)){
        mIIS <- cbind(mIIS[,where.indicators])
      }else{
        stop("'which.ones' not in index")
      }
    }
  } #end if(NROW(x)==1)else..
  return(mIIS)
}
