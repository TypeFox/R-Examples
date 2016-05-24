genTDCM <- function(n, dist, corr, dist.par, model.cens, cens.par, beta, lambda) {
	if (n <= 0) stop("Argument 'n' must be greater than 0", call.=FALSE)
	if ( !( dist %in% c("weibull", "exponential") ) ) stop("Argument 'dist' must be one of 'weibull' or 'exponential'", call.=FALSE)
	if (dist == "weibull") {
		if (corr <= 0 | corr > 1) stop("Argument 'corr' with 'dist=weibull' must be greater than 0 and lower or equal to 1", call.=FALSE)
		if (length(dist.par) != 4) stop("Argument 'dist.par' with 'dist=weibull' must be a vector with lenght 4", call.=FALSE)
		if (dist.par[1] <= 0 | dist.par[2] <= 0 | dist.par[3] <= 0 | dist.par[4] <= 0) stop("Argument 'dist.par' must be greater than 0", call.=FALSE)
	} else if (dist == "exponential") {
		if (corr < -1 | corr > 1) stop("Argument 'corr' with dist='exponential' must be greater or equal to -1 and lower or equal to 1", call.=FALSE)
		if (length(dist.par) != 2) stop("Argument 'dist.par' with 'dist=exponential' must be a vector with lenght 2", call.=FALSE)
		if (dist.par[1] <= 0 | dist.par[2] <= 0) stop("Argument 'dist.par' must be greater than 0", call.=FALSE)
	}
	if ( !( model.cens %in% c("uniform", "exponential") ) ) stop("Argument 'model.cens' must be one of 'uniform' or 'exponential'", call.=FALSE)
	if (cens.par <= 0) stop("Argument 'cens.par' must be greater than 0", call.=FALSE)
	if (length(beta) != 2) stop("Argument 'beta' length must be a vector with length 2", call.=FALSE)
	if (lambda <= 0) stop("Argument 'lambda' must be greater than 0", call.=FALSE)
	b <- dgBIV(n, dist, corr, dist.par)
	if (model.cens == "uniform") {
		rfunc <- runifcens
	} else if (model.cens == "exponential") rfunc <- rexpcens
	mat <- matrix(ncol=6,nrow=1)
	for (k in 1:n) {
		status <- 1
		u <- runif(1, 0, 1)
		z1 <- b[k,2]
		c <- rfunc(1, cens.par)
		if ( u < 1-exp( -lambda*b[k,1]*exp(beta[1]*z1) ) ) {
			t <- -log(1-u)/( lambda*exp(beta[1]*z1) )
			z2 <- 0
		} else {
			t <- -( log(1-u)+lambda*b[k,1]*exp(beta[1]*z1)*( 1-exp(beta[2]) ) )/( lambda*exp(beta[1]*z1+beta[2]) )
			x12 <- b[k,1]
			z2 <- 1
		}
		time <- min(t, c)
		ifelse(t > c, status <- 0, status <- 1)
		if ( u < 1-exp( -lambda*b[k,1]*exp(beta[1]*z1) ) ) {
			aux1 <- c(k, 0, time, status, z1, 0)
			mat <- rbind(mat, aux1)
		} else {
			if (c > x12) {
				aux1 <- c(k, 0, x12, 0, z1, 0)
				mat <- rbind(mat, aux1)
				aux2 <- c(k, x12, time, status, z1, 1)
				mat <- rbind(mat, aux2)
			} else {
				aux1 <- c(k, 0, time, status, z1, 0)
				mat <- rbind(mat, aux1)
			}
		}
	}
	data <- data.frame(mat, row.names=NULL)
	names(data) <- c("id", "start", "stop", "event", "covariate", "tdcov")
	data <- data[-1,]
	row.names(data) <- as.integer( 1:nrow(data) )
	class(data) <- c(class(data), "TDCM")
	return(data)
}
