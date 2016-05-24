genTHMM <- function(n, model.cens, cens.par, beta, covar, rate) {
	if (n <= 0) stop("Argument 'n' must be greater than 0", call.=FALSE)
	if ( !( model.cens %in% c("uniform", "exponential") ) ) stop("Argument 'model.cens' must be one of 'uniform' or 'exponential'", call.=FALSE)
	if (cens.par <= 0) stop("Argument 'cens.par' must be greater than 0", call.=FALSE)
	if (length(beta) != 3) stop("Argument 'beta' length must be a vector with length 3", call.=FALSE)
	if (covar <= 0) stop("Argument 'covar' must be greater than 0", call.=FALSE)
	if (length(rate) != 3) stop("Argument 'rate' must be a vector with length 3", call.=FALSE)
	if (model.cens == "uniform") {
		rfunc <- runifcens
	} else if (model.cens == "exponential") rfunc <- rexpcens
	mat <- matrix(ncol=4, nrow=1)
	for (k in 1:n) {
		z1 <- runif(1, 0, covar)
		c <- rfunc(1, cens.par)
		rate12 <- rate[1]*exp(beta[1]*z1)
		rate13 <- rate[2]*exp(beta[2]*z1)
		rate23 <- rate[3]*exp(beta[3]*z1)
		t12 <- rexp(1, rate12)
		t13 <- rexp(1, rate13)
		t23 <- rexp(1, rate23)
		if (c < min(t12, t13) ) {
			aux1<-c(k,0,1,z1)
			mat<-rbind(mat,aux1)
			aux2<-c(k,c,1,z1)
			mat<-rbind(mat,aux2)
		} else {
			if (t13 < t12) {
				aux1 <- c(k, 0, 1, z1)
				mat <- rbind(mat, aux1)
				aux2 <- c(k, t13, 3, z1)
				mat <- rbind(mat, aux2)
			} else {
				aux1 <- c(k, 0, 1, z1)
				mat <- rbind(mat, aux1)
				aux2 <- c(k, t12, 2, z1)
				mat <- rbind(mat, aux2)
				if (c < t12+t23) {
					aux1 <- c(k, c, 2, z1)
					mat <- rbind(mat, aux1)
				} else {
					aux1 <- c(k ,t12+t23, 3, z1)
					mat <- rbind(mat, aux1)
				}
			}
		}
	}
	data <- data.frame(mat, row.names=NULL)
	names(data) <- c("PTNUM", "time", "state", "covariate")
	data <- data[-1,]
	row.names(data) <- as.integer( 1:nrow(data) )
	class(data) <- c(class(data), "THMM")
	return(data)
}
