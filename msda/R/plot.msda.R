######################################################################
## This function is adapted/modified based on the plot
#   function from
## the glmnet package:
## Jerome Friedman, Trevor Hastie, Robert Tibshirani
#   (2010).
## Regularization Paths for Generalized Linear Models via
#   Coordinate Descent.
##        Journal of Statistical Software, 33(1), 1-22.
##        URL http://www.jstatsoft.org/v33/i01/.

plot.msda <- function(x, xvar = c("norm", "lambda"), ...) {
	theta_nrow = ncol(x$x)
	theta_ncol = x$dim[2]
	nlam <- length(x$lambda)
	ind <- matrix(NA, nlam, theta_ncol)
	tmp <- do.call(cbind,lapply(x$theta,matrix,nrow=theta_nrow,byrow=FALSE))
	par(mfrow = c(1,theta_ncol))
	for(i in seq(theta_ncol)){
		ind <- seq(i,theta_ncol * nlam,by = theta_ncol)
		theta <- tmp[,ind]
	    lambda <- x$lambda
	    df <- x$df
	    xvar <- match.arg(xvar)
	    ##theta should be in 'dgCMatrix' format
	    # which <- nonzero(theta)
	    # theta <- as.matrix(theta[which, ])
	    switch(xvar, norm = {
	        index <- apply(abs(theta), 2, sum)
	        iname <- "L1 Norm"
	    }, lambda = {
	        index <- log(lambda)
	        iname <- "Log Lambda"
	    })
	    xlab <- iname
	    ylab <- "Coefficients"
	    matplot(index, t(theta), lty = 1, xlab = xlab, ylab = ylab, 
	            type = "l", pch = 500, col = rainbow(12, 
	            start = 0.7, end = 0.95), ...) 
	    atdf <- pretty(index)
	    prettydf <- trunc(approx(x = index, y = df, xout = atdf, 
	        rule = 2)$y)
	    axis(3, at = atdf, labels = prettydf, cex.axis = 1, tcl = NA)
	}
} 
