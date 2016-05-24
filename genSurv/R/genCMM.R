genCMM <- function(n, model.cens, cens.par, beta, covar, rate) {
	if (n <= 0) stop("Argument 'n' must be greater than 0", call.=FALSE)
	if ( !( model.cens %in% c("uniform", "exponential") ) ) stop("Argument 'model.cens' must be one of 'uniform' or 'exponential'", call.=FALSE)
	if (cens.par <= 0) stop("Argument 'cens.par' must be greater than 0", call.=FALSE)
	if (length(beta) != 3) stop("Argument 'beta' length must be a vector with length 3", call.=FALSE)
	if (covar <= 0) stop("Argument 'covar' must be greater than 0", call.=FALSE)
	if (length(rate) != 6) stop("Argument 'rate' must be a vector with length 6", call.=FALSE)
	if (model.cens == "uniform") {
		rfunc <- runifcens
	} else if (model.cens == "exponential") rfunc <- rexpcens
	mat <- matrix(ncol=6, nrow=1)
	for (k in 1:n) {
		z1 <- runif(1, 0, covar)
		c <- rfunc(1, cens.par)
		u <- runif(1, 0, 1)
		t12 <- ( -log(1-u)/( rate[1]* exp(beta[1]*z1) ) )^(1/ rate[2])
		u <- runif(1, 0, 1)
		t13 <- ( -log(1-u)/ ( rate[3]* exp(beta[2]*z1) ) )^(1/ rate[4])
		u <- runif(1, 0, 1)
		t23 <- ( -log(1-u)/ ( rate[5]* exp(beta[3]*z1) ) )^(1/ rate[6])
		if (c < min(t12, t13) ) {
			aux1 <- c(k, 0, c, 0, z1, 1)
			mat <- rbind(mat, aux1)
			aux2 <- c(k, 0, c, 0, z1, 2)
			mat <- rbind(mat, aux2)
		} else {
			if (t13 < t12) {
				aux1 <- c(k, 0, t13, 1, z1, 1)
				mat <- rbind(mat, aux1)
				aux2 <- c(k, 0, t13, 0, z1, 2)
				mat <- rbind(mat, aux2)
			} else {
				aux1 <- c(k, 0, t12, 0, z1, 1)
				mat <- rbind(mat, aux1)
				aux2 <- c(k, 0, t12, 1, z1, 2)
				mat <- rbind(mat, aux2)
				if (c < t12+t23) {
					aux1 <- c(k, t12, c, 0, z1, 3)
					mat <- rbind(mat, aux1)
				} else {
					aux1 <- c(k, t12, t12+t23, 1, z1, 3)
					mat <- rbind(mat, aux1)
				}
			}
		}
	}
	data <- data.frame(mat, row.names=NULL)
	data <- data[-1,]
	names(data) <- c("id", "start", "stop", "event", "covariate", "trans")
	row.names(data) <- as.integer( 1:nrow(data) )
	class(data) <- c(class(data), "CMM")
	return(data)
}
