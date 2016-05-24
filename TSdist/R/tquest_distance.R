
TquestDistance <- function(x, y, tx=NULL, ty=NULL, tau) {
  
  # If both temporal indices are missing, equal sampling is assumed and both
  # series begin and end in the same timestamp.
  if (is.null(tx) & is.null(ty)) {
    tx <- seq(0, 1, length.out=length(x))
    ty <- seq(0, 1, length.out=length(y))
  }
  
  # If the temporal index of one of the series is missing, an equal starting and ending point is assumed and the series is sampled constantly in that interval.
  if (is.null(tx)) {
    tx <- seq(ty[1], ty[length(ty)], length.out=length(x))
  }
  if (is.null(ty)) {
    ty <- seq(tx[1], tx[length(tx)], length.out=length(y))
  }
  
  if (class(try(tquestInitialCheck(x, y, tx, ty, tau))) == "try-error") {
    return(NA)
  } else {
   
  # The threshold passing time intervals are identified.
  xbig <- which(x > tau)
  xbig<-tx[xbig]
  ybig <- which(y > tau)
  ybig<-ty[ybig]
  
  # If there are none in both series the resulting distance is 0
  if (length(xbig) == 0 & length(ybig) == 0) {
    d <- 0
  return(d)
  }

  # If on series has threshold passing points and the other does not, 
  # the resulting distance is infinite
  if (length(xbig) == 0 & length(ybig) != 0) {
    d <- Inf
    return(d)
  }
  
  # If on series has threshold passing points and the other does not, 
  # the resulting distance is infinite
  if (length(xbig) != 0 & length(ybig) == 0) {
    d <- Inf
    return(d)
  }

  # The the upper and lower bounds of the threshold trespassing intervals 
  # are found. 
  xinterval <- sort(c(xbig[1], xbig[diff(xbig) != 1], 
                      xbig[which((diff(xbig) != 1)) + 1]))
  
  yinterval<-sort(c(ybig[1], ybig[diff(ybig) != 1], 
                    ybig[which((diff(ybig) != 1)) + 1]))

  # If the number of points obtained in the previous step is uneven, 
  # the last point of the series is included.
  if(length(xinterval) %% 2 == 1) {
    xinterval <- c(xinterval,length(x))
  }
  
  if(length(yinterval) %% 2 == 1) {
    yinterval <- c(yinterval,length(y))
  }


  # The intervals are written in a two column format,
  # the first column is the lower bound and the second column the upper one.
  xinterval <- matrix(xinterval, ncol = 2, byrow=TRUE)
  yinterval <- matrix(yinterval, ncol = 2, byrow=TRUE)

  # The distance between the interval sets of x and the intervals sets of y are 
  # calculated.
  aux1 <- apply(xinterval, 1, function(x) 
    sqrt(rowSums(t(x - t(yinterval)) ^ 2)))

  # For each interval in series x, its closest interval in y is found.
  if(dim(yinterval)[1] != 1) {
    d1 <- sum(apply(aux1, 2, min))
  } else {
    d1 <- min(aux1)
  }

  # The distance between the interval sets of y and the intervals sets of x are 
  # calculated.
  aux2 <- apply(yinterval, 1, function(x) 
    sqrt(rowSums(t(x - t(xinterval)) ^ 2)))
  
  # For each interval in series y, its closest interval in x is found.
  if(dim(xinterval)[1] != 1) {
    d2 <- sum(apply(aux2, 2, min))
  } else { 
    d2 <- min(aux2)
  }

  # The final distance value is calculated.
  d <- 1 / dim(xinterval)[1] * d1 + 1 / dim(yinterval)[1] * d2
  return(d)
  }
}


#  This function checks for possible initial errors: 
tquestInitialCheck <- function(x, y, tx, ty, tau) {
  
  if (! is.numeric(x) | ! is.numeric(y)) {
    stop('The series must be numeric', call.=FALSE)
  }
  if (! is.vector(x) | ! is.vector(y)) {
    stop('The series must be univariate vectors', call.=FALSE)
  }
  if (length(x) <= 1 | length(y) <= 1) {
    stop('The series must have a more than one point', call.=FALSE)
  }
  if (! is.numeric(tau)) {
    stop('The threshold must be numeric', call.=FALSE)
  }
  if (any(is.na(x)) | any(is.na(y))) {
    stop('There are missing values in the series', call.=FALSE)
  } 
  if (any(is.na(tx)) | any(is.na(ty))) {
    stop('There are missing values in the time indices', call.=FALSE)
  }
  
  if (length(tx) != length(x) | length(ty) != length(y)) {
    stop('The length of the time indice must be equal to the length of 
         the series', call.=FALSE)
  }
  
  if (any(tx < 0) | any(ty < 0)) {
      stop('The temporal indices must be positive', call.=FALSE)
    }
  if (any(diff(tx) <= 0) | any(diff(ty) <= 0)) {
      stop('The temporal indices must be strictly increasing', call.=FALSE)
    }
}
