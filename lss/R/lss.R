"lss.eres" <- function(e,z,delta,eps,sh=FALSE)
{
	nobs=length(e)
	ord <- order(e)
	ei <- e[ord]
	zi <- z[ord]
	deltai <- delta[ord]
	tie <- c(diff(ei)>eps,1)
	tie1 <- c(1,diff(ei)>eps)
	dummy <- 1:nobs
	repeats <- diff(c(dummy[c(TRUE,diff(ei)>eps)],nobs+1))
	Ni <- rev(cumsum(rev(zi)))
	di=cumsum(zi*deltai)
	di=di[tie>eps]
	di=c(di[1],diff(di))
	ieb <- 1 - di / Ni[tie1>eps]
	Shat <- cumprod(ieb)
	if(sh)
	{
		return(Shat)	
	}
	Shat <- rep(Shat,repeats)
	edif <- c(diff(ei),0)
	ehat <- rev(cumsum(rev(edif * Shat)))
	ehat[Shat<eps] <- 0
	Shat[Shat<eps] <- 1
	ehat <- ehat/Shat + ei
	eres <- ehat
	eres[dummy[ord]] <- ehat
	eres
}

"lss.betag" <-function(x,y,delta,z)
{
	if(is.vector(x))
		row=1
	else
		row=ncol(x)
	col=ncol(z)

		
	betagm<-matrix(0,ncol=col,nrow=row)
	
	ynew<-1000*(length(y))^2
	if(is.vector(x))
	{
		n1<-length(x)
		n2<- 1
	}
	else
	{
		dimnum<-dim(x)
		n1<-dimnum[1]
		n2<-dimnum[2]
	}
	
	yy0<-rep(y,rep(n1,n1))
	delta1<-rep(delta,rep(n1,n1))
	yy1<-rep(y,n1)
	yy2<-delta1*(yy0-yy1)
	
	xx0<-matrix(rep(as.vector(x),rep(n1,n1*n2)),nrow=n1*n1)
	xx1<-t(matrix(rep(as.vector(t(x)),n1),nrow=n2))
	xx2<-xx0-xx1
	
	for(i in 1:col)
	{
		zz=rep(z[,i],rep(n1,n1))*rep(z[,i],n1)
		xxdif<-xx2*zz*delta1
		xnew<-apply(xxdif,2,sum)
		xnew<-rbind(xxdif,-xnew)
		yynew<-c(yy2*zz,ynew)
	
		if(is.na(LETTERS[c(NA,2)][1])) # if running in R
			fit <- rq(yynew ~ xnew - 1, tau = 0.5)
		else
			fit <- l1fit(xnew, yynew, intercept=FALSE)
		betagm[,i] <- fit$coef

	}

	betagm
}


lss<-function(formula, data, subset, trace=FALSE, mcsize=500, maxiter=10, 
		tolerance=0.001, gehanonly=FALSE, cov=FALSE, na.action=na.exclude)
{
	set.seed(1)
	eps <- .Machine$double.eps^(2/3)

	call <- match.call()
	mf <- match.call(expand.dots = FALSE)
	m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
	mf <- mf[c(1, m)]
	mf$drop.unused.levels <- TRUE
	mf[[1]] <- as.name("model.frame")
	mf <- eval(mf, sys.parent())
	Terms <- attr(mf, "terms")
	xvars <- as.character(attr(Terms, "variables"))
	yvar <- attr(Terms, "response")
	if((yvar <- attr(Terms, "response")) > 0)
		xvars <- xvars[ - yvar]
	else xlevels <- NULL
	y <- model.extract(mf, response)
	x <- model.matrix(Terms, mf)
	x <- as.matrix(x)

	if(all(x[, 1] == 1))
		x <- x[, -1]	

		
	if(ncol(as.matrix(y)) != 2)
		stop("Response must be a right-censored survival object!")

	nobs <- nrow(y)
	if(is.vector(x))
		nvar <- 1
	else
		nvar <- ncol(x)
		

	fit <- list(converged = FALSE, gehanonly=gehanonly, cov=cov, mcsize=mcsize)
	class(fit) <- c("lss")
	fit$call <- call
	fit$nobs <- nobs
	fit$censored <- nobs - sum(y[,2])
	fit$cnames <- dimnames(x)[[2]]
	fit$niter <- 0

	z <- matrix(rexp(nobs*mcsize), ncol=mcsize)
	zdummy <- matrix(rep(1,nobs), ncol=1)

	beta <- lss.betag(x, y[,1], y[,2], zdummy)
	betastar <- lss.betag(x, y[,1], y[,2], z)

	fit$betag <- beta
	
	bbar <- apply(betastar, 1, mean)
	tmp <- betastar - bbar
	fit$gehancov <- tmp %*% t(tmp)/(mcsize - 1)
	fit$gehansd <- sqrt(diag(fit$gehancov))
	fit$gehanzvalue <- beta/fit$gehansd
	fit$gehanpvalue <- (1 - pnorm(abs(fit$gehanzvalue))) * 2
	dimnames(fit$gehancov) <- list(fit$cnames,fit$cnames)
	if(gehanonly)
		return(fit)	

	if(trace)
		cat("\nbetag: ", format(beta), "\n\n")

	niter=0
	if(is.vector(x))
		xbar <- mean(x)
	else
		xbar <- apply(x, 2, mean)
	xm <- x - rep(xbar, rep(nobs, nvar))
	xinv <- solve(t(xm) %*% xm)
	xinvstar <- array(1,dim=c(nvar,nvar,mcsize))
	for(i in 1:mcsize)
	{
		xinvstar[,,i] <- solve(t(xm) %*% (xm*z[,i]))
	}
	
	while(niter < maxiter)
	{
		niter <- niter + 1
		betaprev <- beta

		e <- y[,1] - x %*% beta
		eres <- lss.eres(e, zdummy, y[,2], eps)
		yhat <- y[,2] * y[,1] + (1 - y[,2]) * (eres + x %*% beta)
		ybar <- mean(yhat)
		beta <- xinv %*% (t(xm) %*% (yhat - ybar))

		if(trace)
		{
			cat("Iteration: ", niter)
			cat("\n  Beta: ", format(beta), "\n")
		}

		for(i in 1:mcsize)
		{
			if(is.vector(x))
				e <- y[,1] - x * betastar[,i]
			else
				e <- y[,1] - x %*% betastar[,i]
			eres <- lss.eres(e, z[,i], y[,2], eps)
			if(is.vector(x))
				yhat <- y[,2] * y[,1] + (1 - y[,2]) * (eres + x * betastar[,i])
			else
				yhat <- y[,2] * y[,1] + (1 - y[,2]) * (eres + x %*% betastar[,i])
			ybar <- mean(yhat)
			betastar[,i] <- xinvstar[,,i] %*% (t(xm*z[,i]) %*% (yhat - ybar))
		}
		
		bb <- abs(beta)
		bb[bb<0.01] <- 0.01
		mm <- max(abs(beta - betaprev) / bb)
		if(mm < tolerance)
		{
			fit$residue <- y[,1] - x %*% beta
			fit$km.residue <- lss.eres(fit$residue, zdummy, y[,2], eps, sh=TRUE)
			fit$converged <- TRUE
			break
		}
	}

	if(!fit$converged)
		cat("\nNot converged in", maxiter, "steps\n")
	else 
		cat("\n Converged. Criteria Satisfied: ", tolerance)
	fit$niter <- niter
	fit$lse <- beta
	bbar <- apply(betastar, 1, mean)
	tmp <- betastar - bbar
	fit$tol <- tolerance
	fit$covmatrix <- tmp %*% t(tmp)/(mcsize - 1)
	fit$sd <- sqrt(diag(fit$covmatrix))
	fit$zvalue <- beta/fit$sd
	fit$pvalue <- (1 - pnorm(abs(fit$zvalue))) * 2
	dimnames(fit$covmatrix) <- list(fit$cnames,fit$cnames)

	fit
}

"print.lss" <-function(x, ...)
{
	cat("\nCall:\n")
	dput(x$call)
	cat("\n")
	cat("Number of Observations:  ",x$nobs)
	cat("\n")
	cat("Number of Events:  ", x$nobs-x$censored)
	cat("\n")
	cat("Number of Censored: ", x$censored)  
	cat("\n")
	if(!x$gehanonly){
		cat("Number of Iterations: ", x$niter) 
		cat("\n")
		}
	
	cat("Resampling Number: ", x$mcsize)
	cat("\n\n")

	coef <- array(x$betag, c(length(x$betag), 4))
	dimnames(coef) <- list(x$cnames, c("Estimate", "Std. Error", "Z value", "Pr(>|Z|)"))
	coef[, 2] <- x$gehansd
	coef[, 3] <- x$gehanzvalue
	coef[, 4] <- x$gehanpvalue
	cat("Gehan Estimator:\n")
	print(coef)
	cat("\n")

	if(x$cov)
	{
		cat("Gehan Covariance Matrix:\n")
		print(x$gehancov)
		cat("\n")
	}

	if(!x$gehanonly)
	{
		coef <- array(x$lse, c(length(x$lse), 4))
		dimnames(coef) <- list(x$cnames, c("Estimate", "Std. Error", "Z value", "Pr(>|Z|)"))
		coef[, 2] <- x$sd
		coef[, 3] <- x$zvalue
		coef[, 4] <- x$pvalue
		
		cat("Least-Squares Estimator:\n")
		print(coef)
		cat("\n")
		if(x$cov)
		{
			cat("LSE Covariance Matrix:\n")
			print(x$covmatrix)
			cat("\n")
		}

	}
}

