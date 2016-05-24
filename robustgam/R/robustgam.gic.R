QM.poisson <- function(y, m, c, X, family, w){
    n <- length(y)
    eta <- family$linkfun(m)
    sqrtVar <- sqrt(family$variance(m))
    abeta <- t(X)%*%(expect.poisson(m,c,sqrtVar)*w/sqrtVar*family$mu.eta(eta))/n
    k1 <- m-c*sqrtVar
    k2 <- m+c*sqrtVar
    expect2 <- m*(ppois(k2-2,m)-ppois(ceiling(k1-3),m)) + (1-2*m)*(ppois(k2-1,m)-ppois(ceiling(k1-2),m)) + m*(ppois(k2,m)-ppois(ceiling(k1-1),m))+c^2*(ppois(ceiling(k1-1),m)+1-ppois(k2,m))
    A <- diag(expect2*w^2/family$variance(m)*(family$mu.eta(eta))^2)
    Q <- t(X)%*%A%*%X/n-(abeta)%*%t(abeta)
    expect.cond <- sqrtVar*(ppois(k2-2,m)-ppois(ceiling(k1-3),m)) + (1/sqrtVar-2*sqrtVar)*(ppois(k2-1,m)-ppois(ceiling(k1-2),m)) + sqrtVar*(ppois(k2,m)-ppois(ceiling(k1-1),m)) + c*( dpois(floor(k2),m) + dpois(ceiling(k1-1),m) )
    B <- diag(expect.cond/sqrtVar*w*(family$mu.eta(eta))^2)
    M <- t(X)%*%B%*%X/n
    return(list(Q=Q, M=M))
}

QM.binomial <- function(y, m, c, X, family, w){
    n <- length(y)
    eta <- family$linkfun(m)
    sqrtVar <- sqrt(family$variance(m))
    abeta <- t(X)%*%(expect.binomial(m,c,sqrtVar)*w/sqrtVar*family$mu.eta(eta))/n
    k1 <- m-c*sqrtVar
    k2 <- m+c*sqrtVar
    expect2 <- (Huber.deriv((1-m)/sqrtVar, c))^2*m+(Huber.deriv(-m/sqrtVar, c))^2*(1-m)
    A <- diag(expect2*w^2/family$variance(m)*(family$mu.eta(eta))^2)
    Q <- t(X)%*%A%*%X/n-(abeta)%*%t(abeta)
    expect.cond <- Huber.deriv((1-m)/sqrtVar, c)*(1-m)/family$variance(m)*m + Huber.deriv(-m/sqrtVar, c)*(-m)/family$variance(m)*(1-m)
    B <- diag(expect.cond/sqrtVar*w*(family$mu.eta(eta))^2)
    M <- t(X)%*%B%*%X/n
    return(list(Q=Q, M=M))
}

GIC <- function(y, m, c, X, family, w, sD, gic.constant){
    n <- length(y)
    if (family$family=="poisson") {
	expect <- expect.poisson
	QM <- QM.poisson(y, m, c, X, family, w)
    }
    else if (family$family=="binomial") {
	expect <- expect.binomial
	QM <- QM.binomial(y, m, c, X, family, w)
    }
    if (family$family=="poisson") {
	m <- 1e-6*(m<(1e-6)) + m*(m>=(1e-6))
    }
    if (family$family=="binomial") {
	m <- 1e-6*(m<(1e-6)) + (1-1e-6)*(m>(1-1e-6)) + m*(m>=(1e-6))*(m<=(1-1e-6))
    }
    int1 <- array(dim=n)
    int2 <- array(dim=n)
    for (i in (1:n)){
	int1[i] <- integrate(function(t, c1, y1){Huber.deriv((y1-t)/sqrt(family$variance(t)),c1)/sqrt(family$variance(t))}, lower=y[i], upper=m[i],c1=c, y1=y[i], subdivisions=500, stop.on.error=FALSE)$value
	int2[i] <- integrate(function(t, c1){expect(t,c1,sqrt(family$variance(t)))/sqrt(family$variance(t))}, lower=y[i], upper=m[i],c1=c, subdivisions=500, stop.on.error=FALSE)$value
    }
    Qm <- sum(int1*w)-sum(int2*w)
    gic.comp1 <- -2*Qm
    temp <- svd(QM$M + 2/n*sD)
    gic.comp2 <- sum(diag(temp$v%*%diag(1/temp$d)%*%t(temp$u)%*%(QM$Q)))
    return(list(gic=gic.comp1+gic.constant*gic.comp2, gic.comp1=gic.comp1, gic.comp2=gic.comp2))
}

robustgam.GIC <- function(X, y, family, p=3, K=30, c=1.345, show.msg=FALSE, count.lim=200, w.count.lim=50, smooth.basis="tp", wx=FALSE, sp.min=1e-7, sp.max=1e-3, len=50, show.msg.2=TRUE, gic.constant=log(length(y))){
	X <- matrix(X,nrow=length(y))
	nxs <- ncol(X)
	if (length(sp.min)==1){
		sp.min <- rep(sp.min,nxs)
	}
	if (length(sp.max)==1){
		sp.max <- rep(sp.max,nxs)
	}
	if (length(len)==1){
		len <- rep(len,nxs)
	}
	sp <- list()
	for (j in (1:nxs)){
		sp[[j]] <- exp(seq(log(sp.min[j]), log(sp.max[j]), len=len[j]))
	}
	fit.list <- array(list(),dim=len[1])
	gic <- array(dim=len)
	gic.comp1 <- gic
	gic.comp2 <- gic
	#w <- rep(1, length(y))
	for (i in (1:prod(len))){
		tempind <- arrayInd(i,len)
		tempsp <- sapply(1:nxs,function(jj,sp,ind){sp[[jj]][ind[jj]]},sp=sp,ind=as.vector(tempind))
		fit.list[[i]] <- robustgam(X=X, y=y, family=family, p=p, K=K, c=c, sp=tempsp, show.msg=show.msg, count.lim=count.lim, w.count.lim=w.count.lim, smooth.basis=smooth.basis, wx=wx)
		temp <- GIC(y, fit.list[[i]]$fitted, c, fit.list[[i]]$B, family, fit.list[[i]]$w, fit.list[[i]]$sD, gic.constant)
		gic[tempind] <- temp$gic
		gic.comp1[tempind] <- temp$gic.comp1
		gic.comp2[tempind] <- temp$gic.comp2
		if (show.msg.2) {cat("#### ", as.vector(tempind), ": sp=",tempsp," finished. ####\n")}
	}
	optim.index <- which.min(gic)
	optim.index2 <- arrayInd(optim.index,len)
	optim.fit <- fit.list[[optim.index]]
	optim.sp <- sapply(1:nxs,function(jj,sp,ind){sp[[jj]][ind[jj]]},sp=sp,ind=as.vector(optim.index2))
	return(list(fitted.values=optim.fit$fitted.values, initial.fitted=optim.fit$initial.fitted, beta=optim.fit$beta, optim.index=optim.index, optim.index2=optim.index2, optim.gic=gic[optim.index2], optim.sp=optim.sp, fit.list=fit.list, gic=gic, gic.comp1=gic.comp1, gic.comp2=gic.comp2, sp=sp, gic.constant=gic.constant, optim.fit=optim.fit))
}

robustgam.GIC.optim <- function(X, y, family, p=3, K=30, c=1.345, show.msg=FALSE, count.lim=200, w.count.lim=50, smooth.basis="tp", wx=FALSE, lsp.initial=log(1e-4), lsp.min=-20, lsp.max=10, gic.constant=log(length(y)), method="L-BFGS-B",optim.control=list(trace=1)){
	X <- matrix(X,nrow=length(y))
	nxs <- ncol(X)
	if (length(lsp.initial)==1){
		lsp.initial <- rep(lsp.initial,nxs)
	}
	if (length(lsp.min)==1){
		lsp.min <- rep(lsp.min,nxs)
	}
	if (length(lsp.max)==1){
		lsp.max <- rep(lsp.max,nxs)
	}
	#w <- rep(1, length(y))
	fnr <- function(lsp, X, y, family, p, K, c, show.msg, count.lim, w.count.lim, smooth.basis, wx, gic.constant){
		fit <- robustgam(X=X, y=y, family=family, p=p, K=K, c=c, sp=exp(lsp), show.msg=show.msg, count.lim=count.lim, w.count.lim=w.count.lim, smooth.basis=smooth.basis, wx=wx)
		return(GIC(y=y, m=fit$fitted.values, c=c, X=fit$B, family=family, w=fit$w, sD=fit$sD, gic.constant=gic.constant)$gic)
	}

	gic.optim <- optim(par=lsp.initial, fn=fnr, X=X, y=y, family=family, p=p, K=K, c=c, show.msg=show.msg, count.lim=count.lim, w.count.lim=w.count.lim, smooth.basis=smooth.basis, wx=wx, gic.constant=gic.constant, method=method, lower=lsp.min, upper=lsp.max,control=optim.control)

	optim.fit <- robustgam(X=X, y=y, family=family, p=p, K=K, c=c, sp=exp(gic.optim$par), show.msg=show.msg, count.lim=count.lim, w.count.lim=w.count.lim, smooth.basis=smooth.basis, wx=wx)

	return(list(fitted.values=optim.fit$fitted, beta=optim.fit$beta, beta.fit=optim.fit$betal.fit, gic=gic.optim$value, sp=exp(gic.optim$par), gic.optim=gic.optim, w=optim.fit$w, gic.constant=gic.constant, optim.fit=optim.fit))
}