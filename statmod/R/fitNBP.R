## 	fitNBP.R

fitNBP <- function(y,group=NULL,lib.size=colSums(y),tol=1e-5,maxit=40,verbose=FALSE)
#	Fit multi-group negative-binomial model to SAGE data
#	with Pearson estimation of common overdispersion
#	Gordon Smyth
#	8 July 2006. Last modified 13 July 2009.
{
#	Argument checking
	y <- as.matrix(y)
	if(is.null(group)) group <- rep(1,ncol(y))
	group <- as.factor(group)
	if(length(group) != ncol(y)) stop("length(group) must agree with ncol(y)")

#	Derived quantities
	ngenes <- nrow(y)
	nlib <- ncol(y)
	ngroups <- length(levels(group))
	res.df <- ncol(y)-ngroups
	ind <- matrix(FALSE,nlib,ngroups)
	for (i in 1:ngroups) ind[,i] <- group==levels(group)[i]
	
#	Starting values
	offset <- matrix(1,ngenes,1) %*% log(lib.size)
	mu <- pmax(y,0.5)
	phi <- 0
	w <- mu
	z <- w*(log(mu)-offset)
	beta <- matrix(0,ngenes,ngroups)
	eta <- offset
	for (i in 1:ngroups) {
		beta[,i] <- rowSums(z[,ind[,i],drop=FALSE])/rowSums(w[,ind[,i],drop=FALSE])
		eta[,ind[,i]] <- eta[,ind[,i]]+beta[,i]
	}
	if(verbose) cat("mean coef",colMeans(beta),"\n")
	mu <- exp(eta)

#	Alternating iterations
	iter <- 0
	repeat{
#		Update phi
		iter <- iter+1
		if(iter > maxit) {
			warning("maxit exceeded")
			break
		}
		e2 <- (y-mu)^2
		dV <- mu*mu

#		Need to ensure phi is converging from below
		inneriter <- 0
		repeat {
			inneriter <- inneriter+1
			if(inneriter > 10) stop("problem with inner iteration")
			V <- mu*(1+phi*mu)
			X2 <- sum(e2/V)/res.df-ngenes
			if(X2 >= 0) {
				low <- phi
				break
			} else {
				if(phi==0) break
				if(inneriter > 4)
					phi <- 0.9*phi 
				else
					phi <- (low+phi)/2
				if(verbose) cat("mean disp",phi,"\n")
			}
		}
		if(X2<0) break
		dX2 <- sum(e2/V/V*dV)/res.df
		step.phi <- X2/pmax(dX2,1e-6)
		phi <- phi+step.phi
		conv.crit <- step.phi/(phi+1)
		if(verbose) cat("Conv criterion",conv.crit,"\n")
		if(conv.crit < tol) break
		
#		Update mu
		w <- mu/(1+phi*mu)
		z <- (y-mu)/V*mu
		eta <- offset
		for (i in 1:ngroups) {
			beta[,i] <- beta[,i]+rowSums(z[,ind[,i],drop=FALSE])/rowSums(w[,ind[,i],drop=FALSE])
			eta[,ind[,i]] <- eta[,ind[,i]]+beta[,i]
		}
		if(verbose) cat("mean coef",colMeans(beta),"\n")
		if(verbose) cat("disp",phi,"\n")
		mu <- exp(eta)
	}
	colnames(beta) <- levels(group)
	dimnames(mu) <- dimnames(y)
	list(coefficients=beta,fitted.values=mu,dispersion=phi)
}

