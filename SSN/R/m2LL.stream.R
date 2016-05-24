m2LL.stream <- function(theta, m2LLdata, X, dist.hydro = NULL,
	a.mat = NULL, b.mat = NULL, net.zero = NULL,
	weight = NULL, x.dat = NULL, y.dat = NULL,
	Del.i = NULL, A.5 = NULL,
	CorModels, use.nugget, use.anisotropy, useTailDownWeight,
	EstMeth = "REML",loglik.environment=loglik.environment,
	REs=NULL, scale, maxrang = NULL)
{
	if((max(theta) > 20) | (min(theta) < -20)) return(1e+32)
	theta1 <- untrans.theta(theta = theta, scale = scale)
	if(!is.null(maxrang)) 
		if (any(theta1[!is.na(maxrang)] > maxrang[!is.na(maxrang)])) return(1e+32)
	z <- m2LLdata
	n <- length(z)
	p <- length(X[1,])
	if(length(CorModels) == 0) V <- diag(rep(theta1, times = n)) 
	else V <- makeCovMat(theta1, dist.hydro = dist.hydro, w.matrix = weight,
		a.mat = a.mat, b.mat = b.mat, net.zero = net.zero,
		  x.row = x.dat, y.row = y.dat,
		  x.col = x.dat, y.col = y.dat,
		CorModels = CorModels, useTailDownWeight = useTailDownWeight,
		use.nugget = use.nugget, use.anisotropy = use.anisotropy, REs)
	V.eigenvalues <- eigen(V, symmetric=TRUE, only.values=TRUE)
	if(any(V.eigenvalues$values <=0)) stop("covariance matrix is not positive definite")
	if(!is.null(Del.i)) V <- Del.i*A.5*t((Del.i*A.5) * V)
	qrV <- try(qr(V), silent = T)
	if(class(qrV) == "try-error") return(1e+32)
	ViX <- try(solve(qrV,X), silent = T)
	if(class(ViX) == "try-error") return(1e+32)
	covbi <- crossprod(X,ViX) ## Computationally more efficient than covbi <- t(X) %*% ViX
	covb <- solve(covbi)
	b.hat <- covb %*% crossprod(ViX,z) ##b.hat <- covb %*% t(ViX) %*% z
	r <- z - X %*% b.hat
	f1 <- t(r) %*% solve(qrV,r) + sum(log(abs(diag(qr.R(qrV)))))

	if(EstMeth == "REML") f1 <- f1 + sum(log(svd(covbi)$d))
	nmult <- (n - p)*log(2*pi)
	if(EstMeth == "ML") nmult <- n*log(2*pi)
	result <- f1 + nmult

	## Add this value of theta and the -loglik to an environment for
	## examination later
	RESULT <- get("RESULT",loglik.environment)
	RESULT <- rbind(RESULT,c(theta1,result))
	assign("RESULT",RESULT,loglik.environment)
  
	result

}

