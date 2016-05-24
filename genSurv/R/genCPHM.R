genCPHM <- function(n, model.cens, cens.par, beta, covar) {
	if (n <= 0) stop("Argument 'n' must be greater than 0", call.=FALSE)
	if ( !( model.cens %in% c("uniform", "exponential") ) ) stop("Argument 'model.cens' must be one of 'uniform' or 'exponential'", call.=FALSE)
	if (cens.par <= 0) stop("Argument 'cens.par' must be greater than 0", call.=FALSE)
	if (covar <= 0) stop("Argument 'covar' must be greater than 0", call.=FALSE)
	if (model.cens == "uniform") {
		rfunc <- runifcens
	} else if (model.cens == "exponential") rfunc <- rexpcens
	data <- matrix(data=0, ncol=3, nrow=n)
	for (k in 1:n) {
		status <- 1
		z <- runif(1, 0, covar)
		c <- rfunc(1, cens.par)
		x <- rexp( 1, exp(beta*z) )
		time <- min(x, c)
		ifelse(x > c, status <- 0, status <- 1)
		data[k,] <- c(time, status, z)
	}
	data <- data.frame(data)
	names(data) <- c("time", "status", "covariate")
	row.names(data) <- as.integer( 1:nrow(data) )
	class(data) <- c(class(data), "CPHM")
	return(data)
}
