qmleL <- function(yuima, t, ...){
	

	times <- time(yuima@data@zoo.data[[1]])
	minT <- as.numeric(times[1])
	maxT <- as.numeric(times[length(times)])
	
	if(missing(t) )
	 t <- mean(c(minT,maxT))

	if(t<minT || t>maxT)
	  yuima.stop("time 't' out of bounds")
	grid <- times[which(times<=t)]
	qmle(subsampling(yuima, grid=grid), ...)
}

qmleR <- function(yuima, t, ...){
	
	
	times <- time(yuima@data@zoo.data[[1]])
	minT <- as.numeric(times[1])
	maxT <- as.numeric(times[length(times)])
	
	if(missing(t) )
	t <- mean(c(minT,maxT))
	
	if(t<minT || t>maxT)
	 yuima.stop("time 't' out of bounds")
	grid <- times[which(times>=t)]
	qmle(subsampling(yuima, grid=grid), ...)
}





CPointOld <- function(yuima, param1, param2, print=FALSE, plot=FALSE){
	d.size <- yuima@model@equation.number
	n <- length(yuima)[1]

	d.size <- yuima@model@equation.number

	env <- new.env()
	assign("X",  as.matrix(onezoo(yuima)), envir=env)
	assign("deltaX",  matrix(0, n-1, d.size), envir=env)
	for(t in 1:(n-1))
	 env$deltaX[t,] <- env$X[t+1,] - env$X[t,]
	
	assign("h", deltat(yuima@data@zoo.data[[1]]), envir=env)
	assign("time", as.numeric(index(yuima@data@zoo.data[[1]])), envir=env) 

	QL1 <- pminusquasilogl(yuima=yuima, param=param1, print=print, env)
	QL2 <- pminusquasilogl(yuima=yuima, param=param2, print=print, env)

    D <- sapply(2:(n-1), function(x) sum(QL1[1:x]) + sum(QL2[-(1:x)]))
	D <- c(D[1], D, D[length(D)])
	D <- ts(D, start=0, deltat=deltat(yuima@data@zoo.data[[1]]))
	if(plot)
	 plot(D,type="l", main="change point statistics")
	tau.hat <- index(yuima@data@zoo.data[[1]])[which.min(D)]	

	return(list(tau=tau.hat, param1=param1, param2=param2))
}




# partial quasi-likelihood
# returns a vector of conditionational minus-log-likelihood terms
# the whole negative log likelihood is the sum

pminusquasilogl <- function(yuima, param, print=FALSE, env){
	
	diff.par <- yuima@model@parameter@diffusion
	fullcoef <- diff.par
	npar <- length(fullcoef)
	
	nm <- names(param)
    oo <- match(nm, fullcoef)
    if(any(is.na(oo)))
     oo <- oo[-which(is.na(oo))]
    
    if(any(is.na(oo))) 
		yuima.stop("some named arguments in 'param' are not arguments to the supplied yuima model")

    
    param <- param[order(oo)]
    nm <- names(param)
	
	idx.diff <- match(diff.par, nm)
	
	h <- env$h
	
    theta1 <- unlist(param[idx.diff])
	n.theta1 <- length(theta1)
	n.theta <- n.theta1
	
	d.size <- yuima@model@equation.number
#n <- length(yuima)[1]
	n <- dim(env$X)[1]
	
	vec <- env$deltaX 
	
	K <- -0.5*d.size * log( (2*pi*h) )
	
	

	QL <- 0
	pn <- numeric(n-1)
	diff <- diffusion.term(yuima, param, env)
	dimB <- dim(diff[, , 1])
	
	if(is.null(dimB)){  # one dimensional X
		for(t in 1:(n-1)){
			yB <- diff[, , t]^2
			logdet <- log(yB)
			pn[t] <- K - 0.5*logdet-0.5*vec[t, ]^2/(h*yB) 
			QL <- QL+pn[t]
			
		}
	} else {  # multidimensional X
		for(t in 1:(n-1)){
			yB <- diff[, , t] %*% t(diff[, , t])
			logdet <- log(det(yB))
			if(is.infinite(logdet) ){ # should we return 1e10?
				pn[t] <- log(1)
				yuima.warn("singular diffusion matrix")
				return(1e10)
			}else{
				pn[t] <- K - 0.5*logdet + 
				((-1/(2*h))*t(vec[t, ])%*%solve(yB)%*%vec[t, ]) 
				QL <- QL+pn[t]
			}
		}
	}
	
	if(!is.finite(QL)){
		yuima.warn("quasi likelihood is too small to calculate.")
		QL <- 1e10
	}

	if(print==TRUE){
		yuima.warn(sprintf("NEG-QL: %f, %s", -QL, paste(names(param),param,sep="=",collapse=", ")))
	}
	
	
	return(-pn)
	
}


pminusquasiloglL <- function(yuima, param, print=FALSE, env){
	
	diff.par <- yuima@model@parameter@diffusion
	fullcoef <- diff.par
	npar <- length(fullcoef)
	
	nm <- names(param)
    oo <- match(nm, fullcoef)
    if(any(is.na(oo))) 
	yuima.stop("some named arguments in 'param' are not arguments to the supplied yuima model")
	
    
    param <- param[order(oo)]
    nm <- names(param)
	
	idx.diff <- match(diff.par, nm)
	
	h <- env$h
	
    theta1 <- unlist(param[idx.diff])
	n.theta1 <- length(theta1)
	n.theta <- n.theta1
	
	d.size <- yuima@model@equation.number
	n <- length(yuima)[1]
	
	vec <- env$deltaX 
	
	K <- -0.5*d.size * log( (2*pi*h) )
	
	
	
	QL <- 0
	pn <- numeric(n-1)
	diff <- diffusion.term(yuima, param, env)
	dimB <- dim(diff[, , 1])
	
	if(is.null(dimB)){  # one dimensional X
		for(t in 1:(n-1)){
			yB <- diff[, , t]^2
			logdet <- log(yB)
			pn[t] <- K - 0.5*logdet-0.5*vec[t, ]^2/(h*yB) 
			QL <- QL+pn[t]
			
		}
	} else {  # multidimensional X
		for(t in 1:(n-1)){
			yB <- diff[, , t] %*% t(diff[, , t])
			logdet <- log(det(yB))
			if(is.infinite(logdet) ){ # should we return 1e10?
				pn[t] <- log(1)
				yuima.warn("singular diffusion matrix")
				return(1e10)
			}else{
				pn[t] <- K - 0.5*logdet + 
				((-1/(2*h))*t(vec[t, ])%*%solve(yB)%*%vec[t, ]) 
				QL <- QL+pn[t]
			}
		}
	}
	
	if(!is.finite(QL)){
		yuima.warn("quasi likelihood is too small to calculate.")
		QL <- 1e10
	}
	
	if(print==TRUE){
		yuima.warn(sprintf("NEG-QL: %f, %s", -QL, paste(names(param),param,sep="=",collapse=", ")))
	}
	
	return(-pn)
	
}

# symmetrized version

pminusquasiloglsym <- function(yuima, param, print=FALSE, env){
	
	diff.par <- yuima@model@parameter@diffusion
	fullcoef <- diff.par
	npar <- length(fullcoef)
	
	nm <- names(param)
    oo <- match(nm, fullcoef)
    if(any(is.na(oo))) 
	yuima.stop("some named arguments in 'param' are not arguments to the supplied yuima model")
	
    
    param <- param[order(oo)]
    nm <- names(param)
	
	idx.diff <- match(diff.par, nm)
	
	h <- env$h
	
    theta1 <- unlist(param[idx.diff])
	n.theta1 <- length(theta1)
	n.theta <- n.theta1
	
	d.size <- yuima@model@equation.number

	n <- dim(env$X)[1]-1

    idx0 <- 1:round((n-1)/2)

	vec <- matrix((env$X[2*idx0+1]-2* env$X[2*idx0]+env$X[2*idx0-1])/sqrt(2), length(idx0), dim(env$X)[2])

	K <- -0.5*d.size * log( (2*pi*h) )
	

	QL <- 0
	pn <- numeric(length(idx0))
	diff <- diffusion.term(yuima, param, env)
	dimB <- dim(diff[, , 1])
	
	if(is.null(dimB)){  # one dimensional X
		for(t in idx0){
			yB <- diff[, , 2*t-1]^2
			logdet <- log(yB)
			pn[t] <- K - 0.5*logdet-0.5*vec[t, ]^2/(h*yB) 
			QL <- QL+pn[t]
			
		}
	} else {  # multidimensional X
		for(t in idx0[-length(idx0)]){
			yB <- diff[, , 2*t-1] %*% t(diff[, , 2*t-1])
			logdet <- log(det(yB))
			if(is.infinite(logdet) ){ # should we return 1e10?
				pn[t] <- log(1)
				yuima.warn("singular diffusion matrix")
				return(1e10)
			}else{
				pn[t] <- K - 0.5*logdet + 
				((-1/(2*h))*t(vec[t, ])%*%solve(yB)%*%vec[t, ]) 
				QL <- QL+pn[t]
			}
		}
	}
	
	if(!is.finite(QL)){
		yuima.warn("quasi likelihood is too small to calculate.")
		QL <- 1e10
	}
	
	if(print==TRUE){
		yuima.warn(sprintf("NEG-QL: %f, %s", -QL, paste(names(param),param,sep="=",collapse=", ")))
	}
	

	return(-pn)
	
}


CPoint <- function(yuima, param1, param2, print=FALSE, symmetrized=FALSE, plot=FALSE){
	d.size <- yuima@model@equation.number
	n <- length(yuima)[1]
	
	d.size <- yuima@model@equation.number
	
	env <- new.env()
	assign("X",  as.matrix(onezoo(yuima)), envir=env)
	assign("deltaX",  matrix(0, n-1, d.size), envir=env)
	for(t in 1:(n-1))
	env$deltaX[t,] <- env$X[t+1,] - env$X[t,]
	
	assign("h", deltat(yuima@data@zoo.data[[1]]), envir=env)
	assign("time", as.numeric(index(yuima@data@zoo.data[[1]])), envir=env) 
	
	QL1 <- NULL
	QL2 <- NULL
	
	if(!symmetrized){
	 QL1 <- pminusquasilogl(yuima=yuima, param=param1, print=print, env)
	 QL2 <- pminusquasilogl(yuima=yuima, param=param2, print=print, env)
	} else {
	 QL1 <- pminusquasiloglsym(yuima=yuima, param=param1, print=print, env)
	 QL2 <- pminusquasiloglsym(yuima=yuima, param=param2, print=print, env)
	}
	nn <- length(QL1)
    D <- sapply(2:(nn-1), function(x) sum(QL1[1:x]) + sum(QL2[-(1:x)]))
	D <- c(D[1], D, D[length(D)])
    if(symmetrized)
	 D <- ts(D, start=0, deltat=2*deltat(yuima@data@zoo.data[[1]]))
	else
	 D <- ts(D, start=0, deltat=deltat(yuima@data@zoo.data[[1]]))
	if(plot)
	plot(D,type="l", main="change point statistics")
#	tau.hat <- index(yuima@data@zoo.data[[1]])[which.min(D)]	
	tau.hat <- index(D)[which.min(D)]	
	
	return(list(tau=tau.hat, param1=param1, param2=param2))
}
