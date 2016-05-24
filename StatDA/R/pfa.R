"pfa" <-
function (x, factors, data = NULL, covmat = NULL, n.obs = NA, 
    subset, na.action, start = NULL, scores = c("none", "regression", 
        "Bartlett"), rotation = "varimax", maxiter=5,control = NULL, ...) 
{
#
# Function for Principal Factor Analysis (PFA), P. Filzmoser, 30 March 2004
#
# Input data x should be (robustly) scaled!!!
# covmat (if provided) should be a (robust) covariance OR correlation matrix!!!
#
    sortLoadings <- function(Lambda) {
        cn <- colnames(Lambda)
        Phi <- attr(Lambda, "covariance")
        ssq <- apply(Lambda, 2, function(x) -sum(x^2))
        Lambda <- Lambda[, order(ssq), drop = FALSE]
        colnames(Lambda) <- cn
        neg <- colSums(Lambda) < 0
        Lambda[, neg] <- -Lambda[, neg]
        if (!is.null(Phi)) {
            unit <- ifelse(neg, -1, 1)
            attr(Lambda, "covariance") <- unit %*% Phi[order(ssq), 
                order(ssq)] %*% unit
        }
        Lambda
    }
    cl <- match.call()
    na.act <- NULL
    if (is.list(covmat)) {
        if (any(is.na(match(c("cov", "n.obs"), names(covmat))))) 
            stop("covmat is not a valid covariance list")
        cv <- covmat$cov
        n.obs <- covmat$n.obs
        have.x <- FALSE
    }
    else if (is.matrix(covmat)) {
        cv <- covmat
        have.x <- FALSE
    }
    else if (is.null(covmat)) {
        if (missing(x)) 
            stop("neither x nor covmat supplied")
        have.x <- TRUE
        if (inherits(x, "formula")) {
            mt <- terms(x, data = data)
            if (attr(mt, "response") > 0) 
                stop("response not allowed in formula")
            attr(mt, "intercept") <- 0
            mf <- match.call(expand.dots = FALSE)
            names(mf)[names(mf) == "x"] <- "formula"
            mf$factors <- mf$covmat <- mf$scores <- mf$start <- mf$rotation <- mf$control <- mf$... <- NULL
            mf[[1]] <- as.name("model.frame")
            mf <- eval(mf, parent.frame())
            na.act <- attr(mf, "na.action")
            z <- model.matrix(mt, mf)
        }
        else {
            z <- as.matrix(x)
            if (!missing(subset)) 
                z <- z[subset, , drop = FALSE]
        }
        covmat <- cov.wt(z)
        cv <- covmat$cov
        n.obs <- covmat$n.obs
    }
    else stop("covmat is of unknown type")
    scores <- match.arg(scores)
    if (scores != "none" && !have.x) 
      z <- x
#        stop("requested scores without an x matrix")
    sds <- sqrt(diag(cv))
    cv <- cv/(sds %o% sds)
    p <- ncol(cv)
    dof <- 0.5 * ((p - factors)^2 - p - factors)
#    if (dof < 0) 
#        stop(paste(factors, "factors is too many for", p, "variables"))
    cn <- list(nstart = 1, trace = FALSE, lower = 0.005)
    cn[names(control)] <- control
    more <- list(...)[c("nstart", "trace", "lower", "opt", "rotate")]
    if (length(more)) 
        cn[names(more)] <- more
    if (is.null(start)) {
        start <- (1 - 0.5 * factors/p)/diag(solve(cv))
#        if ((ns <- cn$nstart) > 1) 
#            start <- cbind(start, matrix(runif(ns - 1), p, ns - 
#                1, byrow = TRUE))
    }
    start <- as.matrix(start)
    if (nrow(start) != p) 
        stop(paste("start must have", p, "rows"))
    nc <- ncol(start)
    if (nc < 1) 
        stop("no starting values supplied")
#    best <- Inf
#    for (i in 1:nc) {
#        nfit <- factanal.fit.mle(cv, factors, start[, i], max(cn$lower, 
#            0), cn$opt)
#        if (cn$trace) 
#            cat("start", i, "value:", format(nfit$criteria[1]), 
#                "uniqs:", format(as.vector(round(nfit$uniquenesses, 
#                  4))), "\n")
#        if (nfit$converged && nfit$criteria[1] < best) {
#            fit <- nfit
#            best <- fit$criteria[1]
#        }
#    }
#    if (best == Inf) 
#        stop("Unable to optimize from these starting value(s)")
# PFA:
     fit <- factanal.fit.principal(cv, factors, p=p, start=start[, 1],iter.max=maxiter)
#    cv.red <- cv - diag(1-start)
#    eig <- eigen(cv.red)
#    load <- eig$vectors[,1:factors]%*%diag(sqrt(eig$values[1:factors]))
#    fit <- list(loadings=load,uniqueness=1-start,factors=factors,method="pfa",
#                dof=dof,n.obs=n.obs)
    load <- fit$loadings
    if (rotation != "none") {
        rot <- do.call(rotation, c(list(load), cn$rotate))
        load <- if (is.list(rot)) 
            rot$loadings
        else rot
    }
    fit$loadings <- sortLoadings(load)
    class(fit$loadings) <- "loadings"
    fit$na.action <- na.act
#    if (have.x && scores != "none") {
    if (scores != "none") {
        Lambda <- fit$loadings
#        zz <- scale(z, TRUE, TRUE)
        zz <- z
        switch(scores, regression = {
            sc <- as.matrix(zz) %*% solve(cv, Lambda)
            if (!is.null(Phi <- attr(Lambda, "covariance"))) 
                sc <- sc %*% Phi
        }, Bartlett = {
            d <- 1/fit$uniquenesses
            tmp <- t(Lambda * d)
            sc <- t(solve(tmp %*% Lambda, tmp %*% t(zz)))
        })
        rownames(sc) <- rownames(z)
        colnames(sc) <- colnames(Lambda)
        if (!is.null(na.act)) 
            sc <- napredict(na.act, sc)
        fit$scores <- sc
    }
    if (!is.na(n.obs) && dof > 0) {
        fit$STATISTIC <- (n.obs - 1 - (2 * p + 5)/6 - (2 * factors)/3) * 
            fit$criteria["objective"]
        fit$PVAL <- pchisq(fit$STATISTIC, dof, lower.tail = FALSE)
    }
    fit$n.obs <- n.obs
    fit$call <- cl
    fit
}


"factanal.fit.principal"<-
function(cmat, factors, p = ncol(cmat), start = NULL, iter.max=10 , unique.tol
	 = 0.0001)
{
	dof <- 0.5 * ((p - factors)^2 - p - factors)
	if(dof < 0)
		warning("negative degrees of freedom")
	if(any(abs(diag(cmat) - 1) > .Machine$single.eps))
		stop("must have correlation matrix")
	if(length(start)) {
		if(length(start) != p)
			stop("start is the wrong length")
		if(any(start < 0 | start >= 1))
			stop("all values in start must be between 0 and 1")
		oldcomm <- 1 - start
	}
	else {
		diag(cmat) <- NA
		oldcomm <- apply(abs(cmat), 1, max, na.rm = TRUE)
	}
	diag(cmat) <- oldcomm
	if(iter.max < 0)
		stop("bad value for iter.max")
	ones <- rep(1, factors)
	if(iter.max == 0){
		z <- eigen(cmat, symmetric = TRUE)
		kvals <- z$values[1:factors]
		if(any(kvals <= 0))
			stop("impermissible estimate reached")
		Lambda <- z$vectors[, 1:factors, drop = FALSE] * rep(kvals^0.5, 
			rep.int(p, factors))
                newcomm <- oldcomm
		#newcomm <- Lambda^2 %*% ones
		#diag(cmat) <- newcomm
	}
        if(iter.max >0){
	for(i in 1:iter.max) {
		z <- eigen(cmat, symmetric = TRUE)
		kvals <- z$values[1:factors]
		if(any(kvals <= 0))
			stop("impermissible estimate reached")
		Lambda <- z$vectors[, 1:factors, drop = FALSE] * rep(kvals^0.5, 
			rep.int(p, factors))
		newcomm <- Lambda^2 %*% ones
		if(all(abs(newcomm - oldcomm) < unique.tol)) {
			iter.max <- i
			break
		}
		oldcomm <- newcomm
		diag(cmat) <- newcomm
	}
        }
	dn <- dimnames(cmat)[[1]]
	dimnames(Lambda) <- list(dn, paste("Factor", 1:factors, sep = ""))
	diag(cmat) <- 1
	uniq <- 1 - drop(newcomm)
	names(uniq) <- dn
	ans <- list(loadings = Lambda, uniquenesses = uniq, correlation = cmat, 
		criteria = c(iterations = iter.max), factors = factors, dof = 
		dof, method = "principal")
	class(ans) <- "factanal"
	ans
}
