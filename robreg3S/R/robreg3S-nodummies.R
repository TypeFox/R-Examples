.robreg3S_noDummies <- function(y, x, filter=TRUE, alpha=0.20, xi=0.01, ...){
  x <- as.matrix(x)
  p <- ncol(x)
  n <- nrow(x)
  if( is.null(colnames(x)) ) 
    colnames(x) <- paste0("X",1:p)

  ## filter step
  xf <- x
  if( filter )
  	for(j in 1:p) xf[,j] <- .cfilter.iter(x[,j], alpha=alpha, miter=3)
	if( mean(rowSums(is.na(xf)) > 0) < xi ) xf <- x
  ## scatter estimation step
  init <- GSE(cbind(y, xf), ...)
  ## coeff estimation step
	b1 <- c(solve(init@S[-1,-1]) %*% init@S[-1,1,drop=F] ) 
	a1 <- init@mu[1] - sum(init@mu[-1] * b1)
  
	coeff.3S <- c(a1, b1)
	names(coeff.3S) <- c("(Intercept)", colnames(x))
	acov.3S <- .robreg3S.noDummies.acov(y, init@ximp[,-1], coeff.3S, init@S, init@mu)
	se.3S <- sqrt(diag(acov.3S$asv)/n)
	tab.3S <- cbind( coeff.3S, se.3S, coeff.3S/se.3S, pnorm(abs(coeff.3S/se.3S), lower.tail=FALSE)*2)
	colnames(tab.3S) <- c("Coef", "Asym.Std.Err.", "Z", "Pr(>|Z|)")
	ximp <- init@ximp[,-1]
	colnames(ximp) <- colnames(x)
	result <- list(Summary.Table=tab.3S,
	               coef=coeff.3S,
	               acov=acov.3S$asv,
	               resid=c(y - as.matrix(x)%*%b1 - a1),
	               sigma.hat=acov.3S$se,
	               MD=mahalanobis( cbind(y, x), init@mu, init@S),
	               weight=init@weights,
	               Syx=init@S,
	               myx=init@mu,
	               xfilter=xf,
	               ximpute=ximp
	                 )
	result
}




