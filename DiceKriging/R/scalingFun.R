
##===============================================================
## 'scalingFun1d' computes the 'scaling' transformation at
## the one-dimensional design points x.
##
## The knots and (positive) derivatives definitng the transfo-
## -mation are given in numeric vectors 'knots' and 'eta'.
##
##===============================================================

scalingFun1d <- function(x, knots, eta){
    
  n <- length(x)
  nKnots <- length(knots)
  
  if ( any(x < knots[1]) | any(x > knots[nKnots]) )
    stop("'x' values must be inside the knots")
  
  xs <- sort(x, index.return = TRUE)
  xsorted <- xs$x
  ind <- xs$ix
  
  ## BUG fix for versions > 1.2. 'rightmost.closed'
  ## must be set to TRUE
  inter <- findInterval(x = xsorted, vec = knots,
                        rightmost.closed = TRUE)
  
  inter <- factor(inter, levels = 1:(nKnots-1))
  
  ## iCuts is of length nKnot
  ## the first and last element of 'iCuts' are 0 and  n
  iCuts <- as.numeric(tapply(xsorted, inter, length))
  iCuts[is.na(iCuts)] <- 0
  iCuts <- c(0, cumsum(iCuts))
  scale <- rep(0, n)
  
  ## cat("Calling C\n")
  res <- .C("Scale",
            n = as.integer(n),
            nKnots = as.integer(nKnots),
            iCuts = as.integer(iCuts),  
            x = as.double(xsorted), 
            knots = as.double(knots),
            eta = as.double(eta),
            scale = as.double(scale))
  
  ## CAUTION here: the values are for SORTED x values
  ## they must be put back in the original order.
  newscale <- rep(0, n) 
  newscale[ind] <- res$scale
  
  return(newscale)    
  
}

##==================================================================
## 'scalingFun' applies d one-dimensional 'scaling' transformations
## to n d-dimensional design points
##
## o 'X' must be a n*d matrix,
##
## o 'knots' and 'eta' are list of lentgh d. In both cases,
##    the element i contains the knots/derivatives for the
##    dimension i  
##
## At this time, no control on the list/matrix dimension or length
## 
##==================================================================

scalingFun <- function(X, knots, eta, plot = FALSE) {
 
  d <- NCOL(X)

  transX <- matrix(NA, nrow = NROW(X), ncol = NCOL(X))
  dimnames(transX) <- dimnames(X)
  X <- as.matrix(X)
  
  for (i in 1:d){
    transX[ , i] <- scalingFun1d(x = X[ , i], knots = knots[[i]], eta = eta[[i]])
    #possible bug to correct: transX[ , i] <- scalingFun1d(x = X[ , i], 
    # knots = knots[[colnames(X)[i]]], eta = eta[[colnames(X)[i]]]) 
  }

  if (plot) {
    for(i in 1:d) {
      plot(X[ ,i], transX[ ,i], type = "o",
           pch = 21, col = "orangered", cex = 0.8,
           xlab = "x", ylab = "x 'scaled'")
      abline(v = knots[[i]], col = "SpringGreen3")
    }

  }
  
  return(transX)
  
}


affineScalingFun <- function(X, knots, eta) {
# Here X is meant to be a n*d matrix,
# knots is a (K+1) vector (in this special case),
# and eta is a d*(K+1) matrix of coefficients.
	
	d <- dim(X)[2]

	transformed_X <- matrix(NA, nrow=nrow(X), ncol=ncol(X))

	for(i in seq(1,d)){
		transformed_X[,i] <- scalingFun1d(x=X[,i], knots=knots, eta = eta[i,])
	}

	return(transformed_X)

}
