tim <-
function(x, which.ones=NULL, log.trend=FALSE)
{
  if(NROW(x)==1){
    x.is.scalar <- TRUE
    n <- x
    if(is.null(which.ones)){
      where.indicators <- 2:n
    }else{
      where.indicators <- which(1:n %in% which.ones)
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
    }
  }
  n.where.indicators <- length(where.indicators)
  mTIS <-matrix(0,n,n.where.indicators)
  v1n <- seq(1,n)
  loop.indx <- 1:n.where.indicators
  tmp <- function(i){
    t.trend <- v1n[1:c(n-where.indicators[i]+1)]
    if(log.trend) t.trend <- log(t.trend)
    mTIS[c(where.indicators[i]:n),i] <<- t.trend
  }
  tmp <- sapply(loop.indx,tmp)
  if(x.is.scalar){
    colnames(mTIS) <- paste("tis", where.indicators, sep="")
    mTIS <- as.zoo(mTIS)
  }else{
    colnames(mTIS) <- paste("tis",
      xIndexAsChar[where.indicators], sep="")
    mTIS <- zoo(mTIS, order.by=x.index)
    if(xIsRegular){ mTIS <- as.zooreg(mTIS) }
  }
  return(mTIS)
}
