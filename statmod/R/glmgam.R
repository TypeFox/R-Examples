#  GLMGAM.R

glmgam.fit <- function(X,y,coef.start=NULL,tol=1e-6,maxit=50,trace=FALSE) {
#  Fit gamma generalized linear model with identity link
#  by Levenberg damped Fisher scoring
#  Gordon Smyth
#  12 Mar 2003.  Last revised 3 November 2010.

#  check input
X <- as.matrix(X)
n <- nrow(X)
p <- ncol(X)
if(p > n) stop("More columns than rows in X")
y <- as.vector(y)
if(n != length(y)) stop("length(y) not equal to nrow(X)")
if(n == 0) return(list(coefficients=numeric(0),fitted.values=numeric(0),deviance=numeric(0)))
if(!(all(is.finite(y)) || all(is.finite(X)))) stop("All values must be finite and non-missing")
if(any(y < 0)) stop("y must be non-negative")
maxy <- max(y)
if(maxy==0) return(list(coefficients=rep(0,p),fitted.values=rep(0,n),deviance=NA))
y1 <- pmax(y,maxy*1e-3)

#  starting values
if(is.null(coef.start)) {
	fit <- lm.fit(X,y)
	beta <- fit$coefficients
	mu <- fit$fitted.values
	if(any(mu < 0)) {
		fit <- lm.wfit(X,y,1/y1^2)
		beta <- fit$coefficients
		mu <- fit$fitted.values
	}
	if(any(mu < 0)) {
		fit <- lm.fit(X,rep(mean(y),n))
		beta <- fit$coefficients
		mu <- fit$fitted.values
	}
	if(any(mu < 0)) {
		samesign <- apply(X>0,2,all) | apply(X<0,2,all)
		if(any(samesign)) {
			i <- (1:p)[samesign][1]
			beta <- rep(0,p)
			beta[i] <- lm.wfit(X[,i,drop=FALSE],y,1/y1^2)$coef
			mu <- X[,i] * beta[i]
		} else
			return(list(coefficients=rep(0,p),fitted.values=rep(0,n),deviance=Inf))
	}
} else {
	beta <- coef.start
	mu <- X %*% beta
}
if(any(mu<0)) stop("Starting values give negative fitted values")

deviance.gamma <- function(y,mu) {
	if(any(mu<0)) return(Inf)
	o <- (y < 1e-15) & (mu < 1e-15)
	if(any(o)) {
		if(all(o)) {
			dev <- 0
		} else {
			y1 <- y[!o]
			mu1 <- mu[!o]
			dev <- 2*sum( (y1-mu1)/mu1 - log(y1/mu1) )
		}
	} else {
		dev <- 2*sum( (y-mu)/mu - log(y/mu) )
	}
}

dev <- deviance.gamma(y,mu)

# Scoring iteration with Levenberg damping
iter <- 0
if(trace) cat("Iter =",iter,", Dev =",dev," Beta",beta,"\n")
repeat {
	iter <- iter+1

	# information matrix
	v <- mu^2
	v <- pmax(v,max(v)/10^3)
	XVX <- crossprod(X,vecmat(1/v,X))
	maxinfo <- max(diag(XVX))
	if(iter==1) {
		lambda <- abs(mean(diag(XVX)))/p
		I <- diag(p)
	}

	# score vector
	dl <- crossprod(X,(y-mu)/v)

	# Levenberg damping
	betaold <- beta
	devold <- dev
	lev <- 0
	repeat {
		lev <- lev+1

		# trial step
		R <- chol(XVX + lambda*I)
		dbeta <- backsolve(R,backsolve(R,dl,transpose=TRUE))
		beta <- betaold + dbeta
		mu <- X %*% beta
		dev <- deviance.gamma(y,mu)
		if(dev <= devold || dev/max(mu) < 1e-15) break

		# exit if too much damping
		if(lambda/maxinfo > 1e15) {
			beta <- betaold
			warning("Too much damping - convergence tolerance not achievable")
			break
		}

		# step not successful so increase damping
		lambda <- 2*lambda
		if(trace) cat("Damping increased to",lambda,"\n")
	}

	# iteration output
	if(trace) cat("Iter =",iter,", Dev =",dev," Beta",beta,"\n")

	# keep exiting if too much damping
	if(lambda/maxinfo > 1e15) break

	# decrease damping if successful at first try
	if(lev==1) lambda <- lambda/10

	# test for convergence
	if( crossprod(dl,dbeta) < tol || dev/max(mu) < 1e-15) break

	# test for iteration limit
	if(iter > maxit) break
}

beta <- drop(beta)
names(beta) <- colnames(X)
list(coefficients=beta,fitted.values=as.vector(mu),deviance=dev,iter=iter)
}
