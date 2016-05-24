`cusp.fit` <-
function(y,x=rep(1,length(y)), x.alpha=x, x.beta=x.alpha, weights = rep(1, nobs),
    start=c(rep(0,qr(x.alpha)$rank),rep(1,qr(x.beta)$rank)/qr(x.beta)$rank,c(rep(0,qr(y)$rank-1),1)),
    ..., offset = rep.int(0, nobs), method='L-BFGS-B', optim.method='L-BFGS-B', expected.method = 'delay',
    lower=if(method=='L-BFGS-B') c(rep(-5,qr(x.alpha)$rank+qr(x.beta)$rank),rep(-50,qr(y)$rank),1e-8) else -Inf,
    upper=if(method=='L-BFGS-B') c(rep(+5,qr(x.alpha)$rank+qr(x.beta)$rank),rep(+50,qr(y)$rank), Inf) else +Inf,
    intercept = TRUE){
    x <- as.matrix(x)
    Qr <- qr(x)
    Q  <- scale(zapsmall(qr.Q(Qr)[,1:Qr$rank]),center=FALSE)
    Qscl = attr(Q,'scaled:scale')
    R = diag(c(Qscl,rep(1,NCOL(x)-Qr$rank)), NCOL(x)) %*% qr.R(Qr)
    if(Qr$rank < NCOL(x))
        start <- c(start[1:(2*NCOL(Q))],start[c(length(start)-1,length(start))])
    # regression matrix alpha
    x.alpha <- as.matrix(x.alpha)
    Qr.alpha <- qr(x.alpha)
    ranka <- Qr.alpha$rank; idxa <- 1:ranka;
    Q.alpha  <- scale(zapsmall(qr.Q(Qr.alpha)[,1:ranka]),center=FALSE)
    Qscl.alpha = attr(Q.alpha,'scaled:scale')
    R.alpha = diag(c(Qscl.alpha,rep(1,NCOL(x.alpha)-ranka)), NCOL(x.alpha)) %*% qr.R(Qr.alpha)
    # regression matrix beta
    x.beta <- as.matrix(x.beta)
    Qr.beta <- qr(x.beta)
    rankb <- Qr.beta$rank; idxb <- 1:rankb
    Q.beta  <- scale(zapsmall(qr.Q(Qr.beta)[,1:rankb]),center=FALSE)
    Qscl.beta = attr(Q.beta,'scaled:scale')
    R.beta = diag(c(Qscl.beta,rep(1,NCOL(x.beta)-rankb)), NCOL(x.beta)) %*% qr.R(Qr.beta)
    # dependents matrix
    Y <- as.matrix(y)
    Qr.y <- qr(Y)
    ranky <- Qr.y$rank; idxy <- 1:ranky
    Q.y <- scale(zapsmall(qr.Q(Qr.y)[,1:ranky]),center=FALSE)
    Qscl.y <- attr(Q.y, 'scaled:scale')
    R.y <- diag(c(Qscl.y, rep(1, NCOL(Y) - ranky)), NCOL(Y)) %*% qr.R(Qr.y)
    w <- start[ranka+rankb+1:ranky]
#    w <- c(-w[2:ranky - 1], 1) / w[ranky]
    y <- Q.y %*% backsolve(R.y,w, k = ranky)

    xnames <- dimnames(x)[[2]]
    xnames.alpha <- dimnames(x.alpha)[[2]]
    xnames.beta <- dimnames(x.beta)[[2]]
    xnames.y    <- dimnames(Y)[[2]]
    ynames <- if (is.matrix(y)) rownames(y) else names(y)
    dev.resids <- function (y, mu, wt) if(NCOL(mu)==2) wt * ((cbind(y,y^2/2) - mu)^2) else wt * (y - mu)^2
    nobs <- NROW(y)
    nvars <- NCOL(x.alpha) + NCOL(x.beta) + NCOL(y)
    EMPTY <- nvars == 0
    if (is.null(weights))
        weights <- rep.int(1, nobs)
    if (is.null(offset))
        offset <- rep.int(0, nobs)
    if (EMPTY) {
        eta <- mat.or.vec(nobs, 2) + offset
        colnames(eta) <- c('alpha', 'beta')
        mu <- cusp.expected(eta[,'alpha'], eta[,'beta'], y, method=expected.method)
        dev <- sum(dev.resids(y, mu, weights))
        w <- weights^0.5
        residuals <- if(NCOL(mu)==2) (cbind(y,y^2/2) - mu) else (y - mu)
        residuals <- as.matrix(residuals)
        colnames(residuals) <- colnames(mu)
        good <- rep(TRUE, length(residuals))
        boundary <- conv <- TRUE
        coef <- numeric(0)
        iter <- 0
    }
    else {
        s <- 1

		if(is.loaded("cuspnc")) {
			fit <- optim(start, cusp.nlogLike.c,
    	    	y=Q.y, X.alpha=Q.alpha, X.beta=Q.beta,..., method=method, lower=lower, upper=upper, hessian=TRUE);
		}
		else {
			fit <- optim(start, cusp.nlogLike,
				y=Q.y, X.alpha=Q.alpha, X.beta=Q.beta,..., method=method, lower=lower, upper=upper, hessian=TRUE);
		}
    	w = backsolve(R.y, fit$par[1:ranky+ranka+rankb], k = ranky);
    	sgn = sign(w[ranky]); # sign of a and w are arbitrary, last entry of w should be positive
 # gaat niet goed zo! moet Qr.alpha$pivot en Qr.beta$pivot gebruiken om kolomen uit matrices te kiezen die de R(x.alpha) en R(x.beta) opspannen!!
        ahat <- rep(NA, ncol(x.alpha))
        ahat[Qr.alpha$pivot[idxa]] <- sgn * backsolve(R.alpha, fit$par[idxa], k = ranka);
        bhat <- rep(NA, ncol(x.beta))
        bhat[Qr.beta$pivot[idxb]] <- backsolve(R.beta, fit$par[idxb+ranka], k = rankb);
        what <- rep(NA, ncol(Y))
        what[Qr.y$pivot[idxy]] <- sgn * backsolve(R.y, fit$par[idxy+ranka+rankb], k = ranky);
#    	coef <- coefold <- c( # sign of a and w are arbitrary, last entry of w should be positive
#    	    a = sgn * backsolve(R.alpha, fit$par[1:ranka], k = ranka), rep.int(NA, NCOL(x.alpha)-ranka),
#    	    b = backsolve(R.beta,  fit$par[1:NCOL(Q.beta )+ranka], k = rankb), rep.int(NA, NCOL(x.beta)-rankb),
#    	    w = sgn * backsolve(R.y, fit$par[1:ranky+ranka+rankb], k = ranky), rep.int(NA, NCOL(Y) - ranky))
        coef <- coefold <- c(a = ahat, b = bhat, w = what)
    	RR = matrix(0,ranka+rankb+ranky, ranka+rankb+ranky)
    	RR[1:ranka, 1:ranka] <- R.alpha[1:ranka,1:ranka]
    	RR[1:rankb+ranka, 1:rankb+ranka] <- R.beta[1:rankb,1:rankb]
    	RR[1:ranky + rankb+ranka, 1:ranky + rankb+ranka] <- R.y[1:ranky,1:ranky]
    	fit$Hessian <- t(RR) %*% fit$hessian %*% RR
    	fit <- c(fit, qr(fit$Hessian))
    	fit$rank <- qr(fit$hessian)$rank # Hessian may become quite 'big' which leads to ill rank estimates
    	attr(fit$rank, "ranks") = c(ranka = ranka, rankb = rankb, ranky = ranky)
    	conv <- fit$convergence == 0
#       if (!conv)
#           warning("algorithm did not converge")
    	eps <- 10 * .Machine$double.eps
        xxnames <- c(paste('a[',xnames.alpha[1:ranka],']',sep=''),
                     paste('b[',xnames.beta[1:rankb],']',sep=''),
                     paste('w[',xnames.y[1:ranky],']',sep=''))[fit$pivot]
        eta <- cbind(`alpha`= sgn * Q.alpha %*% fit$par[1:ranka], `beta`=Q.beta %*% fit$par[1:rankb+ranka])
        w <- fit$par[ranka+rankb+1:ranky]
        y <- sgn * Q.y %*% w # the cusp variate that is embeded in the dependents, its sign depends on sgn = sign(w[length(w)])

        colnames(eta) <- c('alpha', 'beta')
        mu <- cusp.expected(eta[,'alpha'], eta[,'beta'], y, method=expected.method)
    	devold <- dev <- s^2 * sum(dev.resids(y, mu, weights))
        residuals <- s * if(NCOL(mu)==2) (cbind(y,y^2/2)-mu) else (y - mu)
        colnames(residuals) <- colnames(mu)
        fit$qr <- as.matrix(fit$qr)
        nr <- nvars - ranka - rankb - ranky
        Rmat <- as.matrix(fit$qr)
        Rmat[row(Rmat) > col(Rmat)] <- 0
        names(coef) <- c(paste('a[',xnames.alpha,']',sep=''), paste('b[',xnames.beta,']',sep=''), paste('w[',xnames.y,']',sep=''))
        colnames(fit$qr) <- xxnames
        dimnames(Rmat) <- list(xxnames, xxnames)
    	fit$start <- start
    	fit$algorithm <- method
    }
    names(residuals) <- ynames
    names(mu) <- ynames
    rownames(eta) <- ynames
    wt <- weights
    names(wt) <- ynames
    names(weights) <- ynames
    names(y) <- ynames
    fit$effects <- double(nobs)
    if (!EMPTY)
        names(fit$effects) <- c(xxnames[seq_len(fit$rank)], rep.int("", nobs - fit$rank))
    wtdmu <- if(intercept) sum(weights * y)/sum(weights) else offset    # (weighted) mean
    nulldev <- s * sum(dev.resids(y, wtdmu, weights)) # (weighted) SS deviances around the (weighted) mean
    n.ok <- nobs - sum(weights == 0)
    nulldf <- n.ok - as.integer(intercept)
    rank <- if (EMPTY) 0 else fit$rank
    resdf <- n.ok - rank
    aic.model <- 2 * fit$value + 2 * rank
	list(coefficients = coef, residuals = residuals, fitted.values = mu,
	    effects = if (!EMPTY) fit$effects, R = if (!EMPTY) Rmat,
	    rank = rank, qr = if (!EMPTY) structure(fit[c("qr", "rank",
	        "qraux", "pivot", "tol")], class = "qr"), family = quasi(),
	    linear.predictors = eta, deviance = dev, aic = aic.model,
	    null.deviance = nulldev, iter = if (!EMPTY) fit$counts[1], weights = wt, prior.weights = weights,
	    df.residual = resdf, df.null = nulldf, y = y, converged = conv,
	    boundary = FALSE, par = fit$par, Hessian = if (!EMPTY) fit$Hessian,
	    hessian.untransformed = if (!EMPTY) fit$hessian,
	    qr.transform= if (!EMPTY) RR, code= if (!EMPTY) fit$convergence, message= if (!EMPTY) fit$message)
}

