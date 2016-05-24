#############################################################
#
#	localdepth.similarity.mahalanobis function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: December, 17, 2008
#	Version: 0.2-5
#
#	Copyright (C) 2008 Claudio Agostinelli and Mario Romanazzi
#
#############################################################

localdepth.similarity.mahalanobis <- function(x, y=NULL, tau, nsamp='all', nmax=1, location=NULL, covariance=NULL, weight=NULL) {
  mahdepth <- function(x, mean, covinv) {
    temp <- (x-mean)
    1/(1+temp%*%covinv%*%temp)
  }
  if (is.null(y))
    y <- x
  if (is.vector(x))
    x <- matrix(x, nrow=1)
  if (is.vector(y))
    y <- matrix(y, nrow=1)
  x <- as.matrix(x)  
  y <- as.matrix(y)
  
  if (ncol(x) < 2) stop('At least two columns must be supplied')
  nx <- nrow(x)
  ny <- nrow(y)
  if (is.null(covariance))
    covariance <- cov(x)
  if (is.null(location))
    location <- as.vector(apply(x, 2, mean))
  if (is.numeric(nsamp) && nsamp <= 0) stop("the argument 'nsamp' must be positive")
  if (is.numeric(nsamp) && nsamp > round(nx*nmax)) {
      warning("Since 'nsamp' is greater than the number of pairwise objects the 'all' method is used")
      nsamp <- 'all'
  }
  if (is.numeric(nsamp)) {
    index <- sample(1:nx, size=nsamp, replace =TRUE)
    nx <- nsamp
    x <- x[index,]
  }

  if (abs(det(covariance)/2) < .Machine$double.eps) stop('The covariance matrix seems to be singular')
  Sinv <- solve(covariance)
  mah <- matrix(0,ny,nx);
  for(i in 1:ny) {
    for(j in 1:nx) { 
      temp <- y[i,]-x[j,]
      mah[i,j] <- sqrt(temp%*%Sinv%*%temp)
    }
  }
  simld <- matrix(0, ny, ny)
  for(i in 1:ny) {
    vetti <- mah[i,]
    for(j in i:ny) {
      vettj <- mah[j,]
      simld[i,j] <- simld[j,i] <- sum(vetti <= tau & vettj <= tau)/nx
    }
  }

  if (!is.null(weight)) {
    mah <- weight(mah)
    simld <- simld*mah
  }
  
  result <- list()
  result$localdepth <- simld
  result$depth <- NA
  result$max.localdepth <- max(result$localdepth)
  result$max.depth <- NA
  result$num <- c(nx,NA)  
  result$call <- match.call()
  result$tau <- tau
  result$x <- x
  result$y <- y
  result$type <- 'exact'
  result$nsamp <- nsamp
  result$method <- 'mahalanobis'
  class(result) <- 'localdepth.similarity'
  return(result)
}
