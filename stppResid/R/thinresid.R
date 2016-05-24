thinresid<-function(X, cifunction = NULL, theta = NULL,  k = NULL, lambda = NULL)
{	
	if(!is.stpp(X))
		stop("X must be an object of type stpp")
	if(is.null(cifunction) && is.null(lambda))
		stop("Either lambda or the conditional intensity function must be specified")
	if(!is.null(lambda)) { 
		if(length(lambda) != length(X$x))
			stop("lambda must be same length as number of points")
		lamb1 <- lambda 
	} else { 
		if(is.null(theta)) {
			lamb1 <- cifunction(X)
		} else 
			lamb1 <- cifunction(X, theta)
	}
	if(is.null(k)) {
		k <- min(lamb1)
		cat("No cutoff rate specified, using minimum lambda_i =", k, "\n")
	}
	if(k <= 0)
		stop("k must be greater than 0")
	ml <- min(lamb1)
	if(k > ml)
		cat("Warning message: \nk is greater than minimum lambda_i = ", ml, ".\nThinned residuals not appropriate.\n" )
	thin.data <- data.frame(cbind(X$x, X$y, X$t))
	prob <- k/lamb1
	u <- runif(length(prob))
	retain <- (u <= prob)
	keep <- thin.data[retain, ]
	deleted <- thin.data[!retain, ]
	names(keep) <- names(deleted) <- c("x", "y", "t")
	y <- list(X = X, k = k, residuals = keep, deleted = deleted)
	class(y) <- "thinresid"
	return(y)
}