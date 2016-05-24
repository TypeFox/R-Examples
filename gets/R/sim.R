sim <-
function(x, which.ones=NULL)
{
  if(NROW(x)==1){
    x.is.scalar <- TRUE
    n <- x
    if(is.null(which.ones)){
      where.indicators <- 2:n
    }else{
      where.indicators <- which.ones
    }
  }else{
    x.is.scalar <- FALSE
    n <- NROW(x)
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
    if(is.null(which.ones)){
      where.indicators <- 2:n
    }else{
      where.indicators <- which(x.index %in% which.ones)
      if(length(where.indicators)==0) stop("'which.ones' not in index")
    } #end if(is.null(which.ones))
  } #end if(NROW(x)==1)else(..)

  n.where.indicators <- length(where.indicators)
  loop.indx <- 1:n.where.indicators
  mSIS <-matrix(0,n,n.where.indicators)
  tmp <- function(i){ mSIS[ c(where.indicators[i]:n) ,i] <<- 1 }
  tmp <- sapply(loop.indx,tmp)
  if(x.is.scalar){
    colnames(mSIS) <- paste("sis", where.indicators, sep="")
    mSIS <- as.zoo(mSIS)
  }else{
    colnames(mSIS) <- paste("sis",
      xIndexAsChar[where.indicators], sep="")
    if(xIsRegular){
      mSIS <- zooreg(mSIS, frequency=xFrequency,
        start=c(xIndexObs[1],xCycle[1]))
    }else{
      mSIS <- zoo(mSIS, order.by=x.index)
    }
  }
  return(mSIS)
}
