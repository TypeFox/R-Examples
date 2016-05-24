hpd <- function(data, alpha = 0.05, ...)
{
    if(!is.data.frame(data))
        data <- read.table(data, header=TRUE)
    nc <- NCOL(data)
	data <- matrix(unlist(data), ncol = nc)
	data <- data[, -1]
    hpd.target <- function(boundary, alpha, fdens, ...){
		#Constraint 1: log f(lower|x) = log f(upper|x) -- operate on log scale (more stable)
		f1 <-  log(fdens(boundary[1],...)) - log(fdens(boundary[2],...))
  
		#Contraint 2: P(param \in C) = 1-alpha (use integrate -- more stable than ppost)
		#f2 <- ppost(boundary[1],...) - ppost(boundary[2])
		#Important: the dpost function has to be vectorized
		f2 <- integrate(fdens,lower=boundary[1],upper=boundary[2],...)$val - (1-alpha)
  
		#Convert multivariate root finding to minimization problem
		#(this can be numerically instable)
		return(f1^2 + f2^2)
	}

	dens <- apply(data, FUN = density, MARGIN = 2)
	fdens <- lapply(dens, FUN = function(obj) {splinefun(obj$x, obj$y)})
	start <- apply(data, MARGIN = 2, FUN = quantile, c(0.5*alpha,1-0.5*alpha))
	opt <- function(j, ...) {
		optim(start[, j], hpd.target,  alpha=alpha, fdens=fdens[[j]], ...)$par
	}
	res <- lapply(1:(nc-1), FUN = opt)
	res <- matrix(unlist(res), ncol = 2, byrow = TRUE)
	res <- data.frame(lower = res[, 1], upper = res[, 2])
	return(res)
}  

hpd.coda <- function(data, alpha = 0.05, ...)
{	
	if(!is.data.frame(data))
        data <- read.table(data, header=TRUE)
	nc <- NCOL(data)
	data <- matrix(unlist(data), ncol = nc)
	data <- data[, -1]
    prob <- 1 - alpha
	data <- coda::as.mcmc(data)
	res <- coda::HPDinterval(data, prob = prob, ...)
  
	return(data.frame(lower = res[,1], upper = res[,2]))
}  
