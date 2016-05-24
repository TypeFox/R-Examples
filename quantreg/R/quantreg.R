
bandwidth.rq <- function(p, n, hs = TRUE, alpha = 0.05)
{
	# Bandwidth selection for sparsity estimation two flavors:
	#	Hall and Sheather(1988, JRSS(B)) rate = O(n^{-1/3})
	#	Bofinger (1975, Aus. J. Stat)  -- rate = O(n^{-1/5})
	# Generally speaking, default method, hs=TRUE is preferred.

	x0 <- qnorm(p)
	f0 <- dnorm(x0)
	if(hs)
            n^(-1/3) * qnorm(1 - alpha/2)^(2/3) *
                ((1.5 * f0^2)/(2 * x0^2 + 1))^(1/3)
	else n^-0.2 * ((4.5 * f0^4)/(2 * x0^2 + 1)^2)^ 0.2
}


plot.rq.process <- function(x, nrow = 3, ncol = 2, ...)
{
    ## Function to plot estimated quantile regression  process
	tdim <- dim(x$sol)
	p <- tdim[1] - 3
	m <- tdim[2]
	oldpar <- par(no.readonly=TRUE)
	par(mfrow = c(nrow, ncol))
	ylab <- dimnames(x$sol)[[1]]
	for(i in 1:p) {
		plot(x$sol[1,], x$sol[3 + i,  ], xlab = "tau",
                     ylab = ylab[3 + i], type = "l")
	}
par(oldpar)
}

print.rqs <- function (x, ...)
{
    if (!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl)
    }
    coef <- coef(x)
    cat("\nCoefficients:\n")
    print(coef, ...)
    rank <- x$rank
    nobs <- nrow(residuals(x))
    p <- nrow(coef)
    rdf <- nobs - p
    cat("\nDegrees of freedom:", nobs, "total;", rdf, "residual\n")
    if (!is.null(attr(x, "na.message")))
        cat(attr(x, "na.message"), "\n")
    invisible(x)
}

"print.rq" <-
function(x, ...)
{
	if(!is.null(cl <- x$call)) {
		cat("Call:\n")
		dput(cl)
	}
	coef <- coef(x)
	cat("\nCoefficients:\n")
	print(coef, ...)
	rank <- x$rank
	nobs <- length(residuals(x))
	if(is.matrix(coef))
		p <- dim(coef)[1]
	else p <- length(coef)
	rdf <- nobs - p
	cat("\nDegrees of freedom:", nobs, "total;", rdf, "residual\n")
	if(!is.null(attr(x, "na.message")))
		cat(attr(x, "na.message"), "\n")
	invisible(x)
}

print.summary.rqs <- function(x, ...) {
    lapply(x, print.summary.rq)
    invisible(x)
}

"print.summary.rq" <-
function(x, digits = max(5, .Options$digits - 2), ...)
{
	cat("\nCall: ")
	dput(x$call)
	coef <- x$coef
	## df <- x$df
	## rdf <- x$rdf
	tau <- x$tau
	cat("\ntau: ")
	print(format(round(tau,digits = digits)), quote = FALSE, ...)
	cat("\nCoefficients:\n")
	print(format(round(coef, digits = digits)), quote = FALSE, ...)
	invisible(x)
}

"rq" <-
function (formula, tau = 0.5, data, subset, weights, na.action, method = "br",
    model = TRUE, contrasts = NULL, ...)
{
    call <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action"), names(mf), 0)
    mf <- mf[c(1,m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval.parent(mf)
    if(method == "model.frame")return(mf)
    mt <- attr(mf, "terms")
    weights <- as.vector(model.weights(mf))
    Y <- model.response(mf)
    X <- model.matrix(mt, mf, contrasts)
    eps <- .Machine$double.eps^(2/3)
    Rho <- function(u,tau) u * (tau - (u < 0))
    if(length(tau)>1){
	if(any(tau < 0) || any(tau > 1))
		stop("invalid tau:  taus should be >= 0 and <= 1")
	if(any(tau == 0)) tau[tau == 0] <- eps
	if(any(tau == 1)) tau[tau == 1] <- 1 -eps
	coef <- matrix(0,ncol(X),length(tau))
        rho <- rep(0,length(tau))
	fitted <- resid <- matrix(0,nrow(X),length(tau))
	for(i in 1:length(tau)){
	    z <- {if (length(weights))
                 	rq.wfit(X, Y, tau = tau[i], weights, method, ...)
                  else rq.fit(X, Y, tau = tau[i], method, ...)
    		 }
	    coef[,i] <- z$coefficients
	    resid[,i] <- z$residuals
            rho[i] <- sum(Rho(z$residuals,tau[i]))
	    fitted[,i] <- Y - z$residuals
	   }
	taulabs <- paste("tau=",format(round(tau,3)))
	dimnames(coef) <- list(dimnames(X)[[2]],taulabs)
	dimnames(resid) <- list(dimnames(X)[[1]],taulabs)
	fit <- z
        fit$coefficients <-  coef
        fit$residuals <- resid
        fit$fitted.values <- fitted
	if(method == "lasso") class(fit) <- c("lassorqs","rqs")
	else if(method == "scad") class(fit) <- c("scadrqs","rqs")
	else class(fit) <- "rqs"
	}
    else{
        process <- (tau < 0 || tau > 1)
	if(tau == 0) tau <- eps
	if(tau == 1) tau <- 1 -eps
        fit <- {
            if (length(weights))
                rq.wfit(X, Y, tau = tau, weights, method, ...)
            else rq.fit(X, Y, tau = tau, method, ...)
           }
	if(process)
		rho <- list(x = fit$sol[1,],y = fit$sol[3,])
	else {
		dimnames(fit$residuals) <- list(dimnames(X)[[1]],NULL)
                rho <-  sum(Rho(fit$residuals,tau))
		}
	if(method == "lasso") class(fit) <- c("lassorq","rq")
	else if(method == "scad") class(fit) <- c("scadrq","rq")
        else class(fit) <- ifelse(process, "rq.process", "rq")
	}
    fit$na.action <- attr(mf, "na.action")
    fit$formula <- formula
    fit$terms <- mt
    fit$xlevels <- .getXlevels(mt,mf)
    fit$call <- call
    fit$tau <- tau
    fit$weights <- weights
    fit$residuals <- drop(fit$residuals)
    fit$rho <- rho
    fit$method <- method
    fit$fitted.values <- drop(fit$fitted.values)

    attr(fit, "na.message") <- attr(m, "na.message")
    if(model) fit$model <- mf
    fit
}
"rq.fit" <-
function(x, y, tau = 0.5, method = "br", ...)
{
    if(length(tau) > 1) {
	    tau <- tau[1]
	    warning("Multiple taus not allowed in rq.fit: solution restricted to first element")
	}

    fit <- switch(method,
		fn = rq.fit.fnb(x, y, tau = tau, ...),
		fnb = rq.fit.fnb(x, y, tau = tau, ...),
		fnc = rq.fit.fnc(x, y, tau = tau, ...),
		pfn = rq.fit.pfn(x, y, tau = tau, ...),
		br = rq.fit.br(x, y, tau = tau, ...),
		lasso = rq.fit.lasso(x, y, tau = tau, ...),
		scad = rq.fit.scad(x, y, tau = tau, ...),
		{
			what <- paste("rq.fit.", method, sep = "")
			if(exists(what, mode = "function"))
				(get(what, mode = "function"))(x, y, ...)
			else stop(paste("unimplemented method:", method))
		}
		)
	fit$fitted.values <- y - fit$residuals
	fit$contrasts <- attr(x, "contrasts")
	fit
}
"rqs.fit"<-
function(x, y, tau = 0.5, tol = 0.0001)
{
# function to compute rq fits for multiple y's
        x <- as.matrix(x)
        p <- ncol(x)
        n <- nrow(x)
        m <- ncol(y)
        z <- .Fortran("rqs",
                as.integer(n),
                as.integer(p),
                as.integer(m),
                as.integer(n + 5),
                as.integer(p + 2),
                as.double(x),
                as.double(y),
                as.double(tau),
                as.double(tol),
                flag = integer(m),
                coef = double(p * m),
                resid = double(n),
                integer(n),
                double((n + 5) * (p + 2)),
                double(n),
                PACKAGE = "quantreg")
        if(sum(z$flag)>0){
                if(any(z$flag)==2)
                        warning(paste(sum(z$flag==2),"out of",m,
                                "BS replications have near singular design"))
                if(any(z$flag)==1)
                        warning(paste(sum(z$flag==1),"out of",m,"may be nonunique"))
                }
        return(t(matrix(z$coef, p, m)))
}
"formula.rq" <-
function (x, ...)
{
    form <- x$formula
    if (!is.null(form)) {
        form <- formula(x$terms)
        environment(form) <- environment(x$formula)
        form
    }
    else formula(x$terms)
}

"predict.rq" <-
function (object, newdata, type = "none", interval = c("none",
    "confidence"), level = 0.95, na.action = na.pass, ...)
{
    if (missing(newdata))
        return(object$fitted)
    else {
        tt <- terms(object)
        Terms <- delete.response(tt)
        m <- model.frame(Terms, newdata, na.action = na.action,
            xlev = object$xlevels)
        if (!is.null(cl <- attr(Terms, "dataClasses")))
            .checkMFClasses(cl, m)
        X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
    }
    pred <- drop(X %*% object$coefficients)
    dots <- list(...)
    if (length(dots$se))
        boot <- (dots$se == "boot")
    else boot <- FALSE
    if (length(dots$mofn))
         mofn <- dots$mofn
    interval <- match.arg(interval)
    if (!interval == "none") {
        if (interval == "confidence") {
            if (type == "percentile") {
                if (boot) {
		    if(exists("mofn")) {# Rescale and recenter!!
		      n <- length(object$fitted)
                      factor <- ifelse(mofn < n, sqrt(mofn/n), 1)
                      XB <- X %*% t(summary(object, cov = TRUE, ...)$B)/factor
                      pl <- apply(XB, 1, function(x) quantile(x, (1 - level)/2))
                      pu <- apply(XB, 1, function(x) quantile(x, 1 - (1 - level)/2))
                      pl <- pred + factor * (pl - pred)
                      pu <- pred + factor * (pu - pred)
                  }
                  else {
                      XB <- X %*% t(summary(object, cov = TRUE, ...)$B)
                      pl <- apply(XB, 1, function(x) quantile(x, (1 - level)/2))
                      pu <- apply(XB, 1, function(x) quantile(x, 1 - (1 - level)/2))
                  }
                  pred <- cbind(pred, pl, pu)
                  colnames(pred) <- c("fit", "lower", "higher")
                }
                else stop("Percentile method requires se = \"boot\".")
            }
            else if (type == "direct") {
                if (boot)
                  stop("Direct method incompatible with bootstrap covariance matrix estimation")
                Z <- rq.fit(object$x, object$y, tau = -1)$sol
                V <- summary(object, cov = TRUE, ...)
                df <- V$rdf
                tfrac <- qt(1 - (1 - level)/2, df)
                Vun <- V$cov * V$scale^2
                tau <- object$tau
                bn <- tfrac * sqrt(diag(X %*% Vun %*% t(X)))
                tauU <- pmin(tau + bn, 1 - 1/df)
                tauL <- pmax(tau - bn, 1/df)
                tauhat <- Z[1, ]
                yhat <- X %*% Z[-(1:3), ]
                n <- nrow(X)
                pl <- yhat[cbind(1:n, cut(tauL, tauhat, label = FALSE))]
                pu <- yhat[cbind(1:n, cut(tauU, tauhat, label = FALSE))]
                pred <- cbind(pred, pl, pu)
                colnames(pred) <- c("fit", "lower", "higher")
            }
            else {
                V <- summary(object, cov = TRUE, ...)
                df <- V$rdf
                tfrac <- qt((1 - level)/2, df)
                sdpred <- sqrt(diag(X %*% V$cov %*% t(X)))
                pred <- cbind(pred, pred + tfrac * sdpred %o%
                  c(1, -1))
                colnames(pred) <- c("fit", "lower", "higher")
            }
        }
        else stop(paste("No interval method for", interval))
    }
    pred
}

"predict.rqs" <-
function (object, newdata, type = "Qhat", stepfun = FALSE, na.action = na.pass, ...)
{
   ## with all defaults
   if(missing(newdata) & !stepfun) return(object$fitted)

   ## otherwise
   tt <- delete.response(terms(object))
   m <- if(missing(newdata)) model.frame(object) else model.frame(tt, newdata,
       na.action = na.action, xlev = object$xlevels)
   if(!is.null(cl <- attr(tt, "dataClasses"))) .checkMFClasses(cl, m)
   X <- model.matrix(tt, m, contrasts = object$contrasts)
   pred <- X %*% object$coefficients

   ## return stepfun or matrix
   if(stepfun) {
       pred <- t(pred)
       taus <- object$tau
       if(type == "Qhat"){
          if(ncol(pred) > 1) {
              f <- as.list(1:ncol(pred))
              for(i in 1:ncol(pred)) f[[i]] <- stepfun(taus, c(pred[1,i], pred[,i]))
          } else
              f <- stepfun(taus, c(pred[1,1], pred[,1]))
          }
       else if(type == "Fhat"){
          if(ncol(pred) > 1) {
              f <- as.list(1:ncol(pred))
              for(i in 1:ncol(pred)) {
                   o <- order(pred[,i])
                   f[[i]] <- stepfun(pred[o,i],c(taus[1],taus[o]))
                   }
          } else
              f <- stepfun(pred[,1],c(taus[1],taus))
          }
       else stop("type must be either 'Qhat' or 'Fhat'")
       return(f)
   } else {
       return(pred)
   }
}


"predict.rq.process" <-
function (object, newdata, type = "Qhat", stepfun = FALSE, na.action = na.pass, ...)
{
    if (missing(newdata))
        return(object$fitted)
    else {
        tt <- terms(object)
        Terms <- delete.response(tt)
        m <- model.frame(Terms, newdata, na.action = na.action,
            xlev = object$xlevels)
        if (!is.null(cl <- attr(Terms, "dataClasses")))
            .checkMFClasses(cl, m)
        X <- model.matrix(Terms, m, contrasts = object$contrasts)
    }
    pred <- t(X %*% object$sol[-(1:3),])
    if(stepfun){
    	taus <- object$sol[1,]
       if(type == "Qhat"){
          if(ncol(pred) > 1) {
              f <- as.list(1:ncol(pred))
              for(i in 1:ncol(pred)) f[[i]] <- stepfun(taus, c(pred[1,i], pred[,i]))
          } else
              f <- stepfun(taus, c(pred[1,1], pred[,1]))
          }
       else if(type == "Fhat"){
          if(ncol(pred) > 1) {
              f <- as.list(1:ncol(pred))
              for(i in 1:ncol(pred)) {
                  o <- order(pred[,i])
                  f[[i]] <- stepfun(pred[,i],c(taus[1],taus))
                  }
          } else
              f <- stepfun(pred[,1],c(taus[1],taus))
          }
       else stop("type must be either 'Qhat' or 'Fhat'")
       return(f)
	}
    return(pred)
}
"rearrange"  <- function (f, xmin, xmax)
# Revised Version September 11 2007.
{
    if (is.list(f))
        lapply(f, rearrange)
    else {
        if (!is.stepfun(f))
            stop("Only stepfuns can be rearranged.\n")
    	call	<- attributes(f)$call;
    	right <- call[match("right",names(call))]=="TRUE()"
        x 	<- knots(f)
	n	<- length(x)
	if(missing(xmin)) xmin = x[1]
	if(missing(xmax)) xmax = x[n]
	x	<- x[(x >= xmin) & (x <= xmax)]
	x	<- c(xmin, x, xmax)
	n	<- length(x)
	y	<- f(x)
        o 	<- ifelse(rep(right,n-1), order(y[-1])+1, order(y[-n]))
        x 	<- cumsum(c(x[1], diff(x)[o - right]))
        y 	<- y[o]
        y 	<- c(y[1], y, max(y))
        stepfun(x, y, right = right)
    }

}



# Function to compute regression quantiles using original simplex approach
# of Barrodale-Roberts/Koenker-d'Orey.  There are several options.
# The options are somewhat different than those available for the Frisch-
# Newton version of the algorithm, reflecting the different natures of the
# problems typically solved.  Succintly BR for "small" problems, FN for
# "large" ones.  Obviously, these terms are conditioned by available hardware.
#
# Basically there are two modes of use:
# 1.  For Single Quantiles:
#
#       if tau is between 0 and 1 then only one quantile solution is computed.
#
#       if ci = FALSE  then just the point estimate and residuals are returned
#		If the column dimension of x is 1 then ci is set to FALSE since
#		since the rank inversion method has no proper null model.
#       if ci = TRUE  then there are two options for confidence intervals:
#
#               1.  if iid = TRUE we get the original version of the rank
#                       inversion intervals as in Koenker (1994)
#               2.  if iid = FALSE we get the new version of the rank inversion
#                       intervals which accounts for heterogeneity across
#                       observations in the conditional density of the response.
#                       The theory of this is described in Koenker-Machado(1999)
#               Both approaches involve solving a parametric linear programming
#               problem, the difference is only in the factor qn which
#               determines how far the PP goes.  In either case one can
#               specify two other options:
#                       1. interp = FALSE returns two intervals an upper and a
#                               lower corresponding to a level slightly
#                               above and slightly below the one specified
#                               by the parameter alpha and dictated by the
#                               essential discreteness in the test statistic.
#				interp = TRUE  returns a single interval based on
#                               linear interpolation of the two intervals
#                               returned:  c.values and p.values which give
#                               the critical values and p.values of the
#                               upper and lower intervals. Default: interp = TRUE.
#                       2.  tcrit = TRUE uses Student t critical values while
#                               tcrit = FALSE uses normal theory ones.
# 2. For Multiple Quantiles:
#
#       if tau < 0 or tau >1 then it is presumed that the user wants to find
#       all of the rq solutions in tau, and the program computes the whole
#	quantile regression solution as a process in tau, the resulting arrays
#	containing the primal and dual solutions, betahat(tau), ahat(tau)
#       are called sol and dsol.  These arrays aren't printed by the default
#       print function but they are available as attributes.
#       It should be emphasized that this form of the solution can be
#	both memory and cpu quite intensive.  On typical machines it is
#	not recommended for problems with n > 10,000.
#	In large problems a grid of solutions is probably sufficient.
#
rq.fit.br <-
function (x, y, tau = 0.5, alpha = 0.1, ci = FALSE, iid = TRUE,
	interp = TRUE, tcrit = TRUE)
{
    tol <- .Machine$double.eps^(2/3)
    eps <- tol
    big <- .Machine$double.xmax
    x <- as.matrix(x)
    p <- ncol(x)
    n <- nrow(x)
    ny <- NCOL(y)
    nsol <- 2
    ndsol <- 2
    # Check for Singularity of X since br fortran isn't very reliable about this
    if (qr(x)$rank < p)
        stop("Singular design matrix")
    if (tau < 0 || tau > 1) {
        nsol <- 3 * n
        ndsol <- 3 * n
        lci1 <- FALSE
        qn <- rep(0, p)
        cutoff <- 0
        tau <- -1
    }
    else {
        if (p == 1)
            ci <- FALSE
        if (ci) {
            lci1 <- TRUE
            if (tcrit)
                cutoff <- qt(1 - alpha/2, n - p)
            else cutoff <- qnorm(1 - alpha/2)
            if (!iid) {
                h <- bandwidth.rq(tau, n, hs = TRUE)
                bhi <- rq.fit.br(x, y, tau + h, ci = FALSE)
                bhi <- coefficients(bhi)
                blo <- rq.fit.br(x, y, tau - h, ci = FALSE)
                blo <- coefficients(blo)
                dyhat <- x %*% (bhi - blo)
                if (any(dyhat <= 0)) {
                  pfis <- (100 * sum(dyhat <= 0))/n
                  warning(paste(pfis, "percent fis <=0"))
                }
                f <- pmax(eps, (2 * h)/(dyhat - eps))
                qn <- rep(0, p)
                for (j in 1:p) {
                  qnj <- lm(x[, j] ~ x[, -j] - 1, weights = f)$resid
                  qn[j] <- sum(qnj * qnj)
                }
            }
            else qn <- 1/diag(solve(crossprod(x)))
        }
        else {
            lci1 <- FALSE
            qn <- rep(0, p)
            cutoff <- 0
        }
    }
    z <- .Fortran("rqbr", as.integer(n), as.integer(p), as.integer(n +
        5), as.integer(p + 3), as.integer(p + 4), as.double(x),
        as.double(y), as.double(tau), as.double(tol), flag = as.integer(1),
        coef = double(p), resid = double(n), integer(n), double((n +
            5) * (p + 4)), double(n), as.integer(nsol), as.integer(ndsol),
        sol = double((p + 3) * nsol), dsol = double(n * ndsol),
        lsol = as.integer(0), h = integer(p * nsol), qn = as.double(qn),
        cutoff = as.double(cutoff), ci = double(4 * p), tnmat = double(4 *
            p), as.double(big), as.logical(lci1), PACKAGE = "quantreg")
    if (z$flag != 0)
        warning(switch(z$flag, "Solution may be nonunique", "Premature end - possible conditioning problem in x"))
    if (tau < 0 || tau > 1) {
        sol <- matrix(z$sol[1:((p + 3) * z$lsol)], p + 3)
        dsol <- matrix(z$dsol[1:(n * z$lsol)], n)
        vnames <- dimnames(x)[[2]]
        dimnames(sol) <- list(c("tau", "Qbar", "Obj.Fun", vnames),
            NULL)
        return(list(sol = sol, dsol = dsol))
    }
    if (!ci) {
        coef <- z$coef
        dual <- z$dsol[1:n]
        names(coef) <- dimnames(x)[[2]]
        return(list(coefficients = coef, x = x, y = y, residuals = y - x %*% z$coef,
		dual = dual))
    }
    if (interp) {
        Tn <- matrix(z$tnmat, nrow = 4)
        Tci <- matrix(z$ci, nrow = 4)
        Tci[3, ] <- Tci[3, ] + (abs(Tci[4, ] - Tci[3, ]) * (cutoff -
            abs(Tn[3, ])))/abs(Tn[4, ] - Tn[3, ])
        Tci[2, ] <- Tci[2, ] - (abs(Tci[1, ] - Tci[2, ]) * (cutoff -
            abs(Tn[2, ])))/abs(Tn[1, ] - Tn[2, ])
        Tci[2, ][is.na(Tci[2, ])] <- -big
        Tci[3, ][is.na(Tci[3, ])] <- big
        coefficients <- cbind(z$coef, t(Tci[2:3, ]))
        vnames <- dimnames(x)[[2]]
        cnames <- c("coefficients", "lower bd", "upper bd")
        dimnames(coefficients) <- list(vnames, cnames)
        residuals <- y - drop(x %*% z$coef)
        return(list(coefficients = coefficients, residuals = residuals))
    }
    else {
        Tci <- matrix(z$ci, nrow = 4)
        coefficients <- cbind(z$coef, t(Tci))
        residuals <- y - drop(x %*% z$coef)
        vnames <- dimnames(x)[[2]]
        cnames <- c("coefficients", "lower bound", "Lower Bound",
            "upper bd", "Upper Bound")
        dimnames(coefficients) <- list(vnames, cnames)
        c.values <- t(matrix(z$tnmat, nrow = 4))
        c.values <- c.values[, 4:1]
        dimnames(c.values) <- list(vnames, cnames[-1])
        p.values <- if (tcrit)
            matrix(pt(c.values, n - p), ncol = 4)
        else matrix(pnorm(c.values), ncol = 4)
        dimnames(p.values) <- list(vnames, cnames[-1])
        list(coefficients = coefficients, residuals = residuals,
            c.values = c.values, p.values = p.values)
    }
}

"rq.fit.fnb" <-
function (x, y, tau = 0.5, beta = 0.99995, eps = 1e-06)
{
    n <- length(y)
    p <- ncol(x)
    if (n != nrow(x))
        stop("x and y don't match n")
    if (tau < eps || tau > 1 - eps)
        stop("No parametric Frisch-Newton method.  Set tau in (0,1)")
    rhs <- (1 - tau) * apply(x, 2, sum)
    d   <- rep(1,n)
    u   <- rep(1,n)
    wn <- rep(0,10*n)
    wn[1:n] <- (1-tau) #initial value of dual solution
    z <- .Fortran("rqfnb", as.integer(n), as.integer(p), a = as.double(t(as.matrix(x))),
        c = as.double(-y), rhs = as.double(rhs), d = as.double(d),as.double(u),
        beta = as.double(beta), eps = as.double(eps),
        wn = as.double(wn), wp = double((p + 3) * p),
        it.count = integer(3), info = integer(1),PACKAGE= "quantreg")
    if (z$info != 0)
        stop(paste("Error info = ", z$info, "in stepy: singular design"))
    coefficients <- -z$wp[1:p]
    names(coefficients) <- dimnames(x)[[2]]
    residuals <- y - x %*% coefficients
    list(coefficients=coefficients, tau=tau, residuals=residuals)
}

"rq.fit.fnc" <-
function (x, y, R, r, tau = 0.5, beta = 0.9995, eps = 1e-06)
{
    n1 <- length(y)
    n2 <- length(r)
    p <- ncol(x)
    if (n1 != nrow(x))
        stop("x and y don't match n1")
    if (n2 != nrow(R))
        stop("R and r don't match n2")
    if (p != ncol(R))
        stop("R and x don't match p")
    if (tau < eps || tau > 1 - eps)
        stop("No parametric Frisch-Newton method.  Set tau in (0,1)")
    rhs <- (1 - tau) * apply(x, 2, sum)
    u <- rep(1, max(n1,n2)) #upper bound vector and scratch vector
    wn1 <- rep(0, 9 * n1)
    wn1[1:n1] <- (1 - tau) #store the values of x1
    wn2 <- rep(0, 6 * n2)
    wn2[1:n2] <- 1 #store the values of x2
    z <- .Fortran("rqfnc", as.integer(n1), as.integer(n2), as.integer(p),
	a1 = as.double(t(as.matrix(x))), c1 = as.double(-y),
	a2 = as.double(t(as.matrix(R))), c2 = as.double(-r),
	rhs = as.double(rhs), d1 = double(n1), d2 = double(n2),
	as.double(u), beta = as.double(beta), eps = as.double(eps),
        wn1 = as.double(wn1), wn2 = as.double(wn2), wp = double((p + 3) * p),
        it.count = integer(3), info = integer(1), PACKAGE = "quantreg")
    if (z$info != 0)
        stop(paste("Error info = ", z$info, "in stepy2: singular design"))
    coefficients <- -z$wp[1:p]
    names(coefficients) <- dimnames(x)[[2]]
    residuals <- y - x %*% coefficients
    it.count <- z$it.count
    list(coefficients=coefficients, tau=tau, residuals=residuals, it = it.count)
}
"rq.fit.scad" <-
function (x, y, tau = 0.5, alpha = 3.2, lambda = 1, start = "rq", beta = 0.9995, eps = 1e-06)
{
    n <- length(y)
    p <- ncol(x)
    if (n != nrow(x))
        stop("x and y don't match n")
    if (tau < eps || tau > 1 - eps)
        stop("No parametric Frisch-Newton method.  Set tau in (0,1)")
    if(length(lambda) == 1)
         lambda <- c(0,rep(lambda,p-1))
    if(length(lambda) != p)
          stop(paste("lambda must be either of length ",p," or length one"))
    if(any(lambda < 0))
          stop("negative lambdas disallowed")
    R <- diag(lambda,nrow = length(lambda))
    R <- R[which(lambda != 0),, drop = FALSE]
    r <- rep(0,nrow(R))
    X <- rbind(x,  R)
    Y <- c(y, r)
    N <- length(Y)
    rhs <- (1 - tau) * apply(x, 2, sum) + apply(R,2,sum)
    dscad <- function(x, a = 3.7, lambda = 2){
        lambda *  sign(x) *  (abs(x) <= lambda) +
        sign(x) * (a * lambda -  abs(x)) / (a - 1) *
                (abs(x) <= a * lambda) * (abs(x) > lambda)
        }
    binit <- switch(start,
    	rq =   rq.fit.fnb(x, y, tau = tau)$coef[-1],
    	lasso = rq.fit.lasso(x, y, tau = tau, lambda = lambda)$coef[-1]
	)
    coef <- rep(.Machine$double.xmax,p)
    vscad <- rhs - c(0,dscad(binit) * sign(binit))
    it <- 0
    while(sum(abs(binit - coef[-1])) > eps){
	it <- it + 1
    	d <- rep(1, N)
    	u <- rep(1, N)
    	wn <- rep(0, 10 * N)
    	wn[1:N] <- c(rep((1 - tau),n),rep(.5,nrow(R)))
	vrhs <- rhs - vscad
	binit <- coef[-1]
    	z <- .Fortran("rqfnb", as.integer(N), as.integer(p), a = as.double(t(as.matrix(X))),
        	c = as.double(-Y), vrhs = as.double(vrhs), d = as.double(d),
        	as.double(u), beta = as.double(beta), eps = as.double(eps),
        	wn = as.double(wn), wp = double((p + 3) * p), 
            	it.count = integer(3), info = integer(1), PACKAGE = "quantreg")
	coef <- -z$wp[1:p]
	vscad <- c(0,dscad(coef[2:p]) * sign(coef[2:p]))
	}
    if (z$info != 0)
        stop(paste("Error info = ", z$info, "in stepy2: singular design"))
    coefficients <- -z$wp[1:p]
    names(coefficients) <- dimnames(x)[[2]]
    residuals <- y - x %*% coefficients
    it.count <- z$it.count
    list(coefficients=coefficients, residuals=residuals, tau = tau,
		lambda = lambda, it = it.count)
}


"rq.fit.lasso" <-
function (x, y, tau = 0.5, lambda = 1, beta = 0.9995, eps = 1e-06)
{
    n <- length(y)
    p <- ncol(x)
    if (n != nrow(x))
        stop("x and y don't match n")
    if(length(lambda) == 1)
         lambda <- c(0,rep(lambda,p-1))
    if(length(lambda) != p)
          stop(paste("lambda must be either of length ",p," or length one"))
    if(any(lambda < 0))
          stop("negative lambdas disallowed")
    R <- diag(lambda,nrow = length(lambda))
    R <- R[which(lambda != 0),, drop = FALSE]
    r <- rep(0,nrow(R))
    if (tau < eps || tau > 1 - eps)
        stop("No parametric Frisch-Newton method.  Set tau in (0,1)")
    X <- rbind(x, R)
    Y <- c(y, r)
    N <- length(Y)
    rhs <- (1 - tau) * apply(x, 2, sum) + 0.5 * apply(R,2,sum)
    d <- rep(1, N)
    u <- rep(1, N)
    wn <- rep(0, 10 * N)
    wn[1:N] <- rep(0.5, N)
    z <- .Fortran("rqfnb", as.integer(N), as.integer(p), a = as.double(t(as.matrix(X))),
        c = as.double(-Y), rhs = as.double(rhs), d = as.double(d),
        as.double(u), beta = as.double(beta), eps = as.double(eps),
        wn = as.double(wn), wp = double((p + 3) * p), aa = double(p *
            p), it.count = integer(3), info = integer(1), PACKAGE = "quantreg")
    if (z$info != 0)
        stop(paste("Error info = ", z$info, "in stepy2: singular design"))
    coefficients <- -z$wp[1:p]
    names(coefficients) <- dimnames(x)[[2]]
    residuals <- y - x %*% coefficients
    it.count <- z$it.count
    list(coefficients=coefficients, residuals=residuals, tau = tau,
	lambda = lambda, it = it.count)
}

"rq.fit.pfn" <-
# This is an implementation (purely in R) of the preprocessing phase
# of the rq algorithm described in Portnoy and Koenker, Statistical
# Science, (1997) 279-300.  In this implementation it can be used
# as an alternative method for rq() by specifying method="pfn"
# It should probably be used only on very large problems and then
# only with some caution.  Very large in this context means roughly
# n > 100,000.  The options are described in the paper, and more
# explicitly in the code.  Again, it would be nice perhaps to have
# this recoded in a lower level language, but in fact this doesn't
# seem to make a huge difference in this case since most of the work
# is already done in the rq.fit.fnb calls.
#
function(x, y, tau = 0.5,  Mm.factor = 0.8,
	max.bad.fixup = 3, eps = 1e-6)
{
	#rq function for n large --
	n <- length(y)
	if(nrow(x) != n)
		stop("x and y don't match n")
	if(tau < 0 | tau > 1)
		stop("tau outside (0,1)")
	p <- ncol(x)
	m <- round(((p + 1) * n)^(2/3))
	not.optimal <- TRUE
	while(not.optimal) {
		if(m < n)
			s <- sample(n, m)
		else {
			b <- rq.fit.fnb(x, y, tau = tau,  eps = eps)$coef
			break
		}
		xx <- x[s,  ]
		yy <- y[s]
		z <- rq.fit.fnb(xx, yy, tau = tau,  eps = eps)
		xxinv <- solve(chol(crossprod(xx)))
		band <- sqrt(((x %*% xxinv)^2) %*% rep(1, p))
		#sqrt(h<-ii)
		r <- y - x %*% z$coef
		M <- Mm.factor * m
		lo.q <- max(1/n, tau - M/(2 * n))
		hi.q <- min(tau + M/(2 * n), (n - 1)/n)
		kappa <- quantile(r/pmax(eps, band), c(lo.q, hi.q))
		sl <- r < band * kappa[1]
		su <- r > band * kappa[2]
		bad.fixup <- 0
		while(not.optimal & (bad.fixup < max.bad.fixup)) {
			xx <- x[!su & !sl,  ]
			yy <- y[!su & !sl]
			if(any(sl)) {
				glob.x <- c(t(x[sl,  , drop = FALSE]) %*% rep(
					1, sum(sl)))
				glob.y <- sum(y[sl])
				xx <- rbind(xx, glob.x)
				yy <- c(yy, glob.y)
			}
			if(any(su)) {
				ghib.x <- c(t(x[su,  , drop = FALSE]) %*% rep(
					1, sum(su)))
				ghib.y <- sum(y[su])
				xx <- rbind(xx, ghib.x)
				yy <- c(yy, ghib.y)
			}
			z <- rq.fit.fnb(xx, yy, tau = tau,  eps = eps)
			b <- z$coef
			r <- y - x %*% b
			su.bad <- (r < 0) & su
			sl.bad <- (r > 0) & sl
			if(any(c(su.bad, sl.bad))) {
				if(sum(su.bad | sl.bad) > 0.10000000000000001 *
					M) {
					warning("Too many fixups:  doubling m")
					m <- 2 * m
					break
				}
				su <- su & !su.bad
				sl <- sl & !sl.bad
				bad.fixup <- bad.fixup + 1
			}
			else not.optimal <- FALSE
		}
	}
	coefficients <- b
	names(coefficients) <- dimnames(x)[[2]]
	residuals <- y - x %*% b
	return(list(coefficients=coefficients, tau=tau,
		residuals=residuals))
}

"rq.wfit" <-
function(x, y, tau = 0.5, weights, method = "br",  ...)
{
	if(any(weights < 0))
		stop("negative weights not allowed")
	if(length(tau) > 1) {
	    tau <- tau[1]
	    warning("Multiple taus not allowed in rq.wfit: solution restricted to first element")
	}
	contr <- attr(x, "contrasts")
	wx <- x * weights
	wy <- y * weights
	fit <- switch(method,
		fn = rq.fit.fnb(wx, wy, tau = tau, ...),
		fnb = rq.fit.fnb(wx, wy, tau = tau, ...),
		br = rq.fit.br(wx, wy, tau = tau, ...),
		fnc = rq.fit.fnc(wx, wy, tau = tau, ...),
                pfn = rq.fit.pfn(wx, wy, tau = tau, ...), {
			what <- paste("rq.fit.", method, sep = "")
			if(exists(what, mode = "function"))
				(get(what, mode = "function"))(x, y, ...)
			else stop(paste("unimplemented method:", method))
		}
		)
        if(length(fit$sol))
            fit$fitted.values <- x %*% fit$sol[-(1:3),]
        else
            fit$fitted.values <- x %*% fit$coef
	fit$residuals <- y - fit$fitted.values
	fit$contrasts <- attr(x, "contrasts")
	fit$weights <- weights
	fit
}

"summary.rqs" <-
function (object, ...) {
        taus <- object$tau
        xsum <- as.list(taus)
        for(i in 1:length(taus)){
                xi <- object
                xi$coefficients <- xi$coefficients[,i]
                xi$residuals <- xi$residuals[,i]
                xi$tau <- xi$tau[i]
                class(xi) <- "rq"
                xsum[[i]] <- summary(xi, ...)
		if(class(object)[1] == "dynrqs"){
		    class(xsum[[1]]) <- c("summary.dynrq", "summary.rq")
	            if(i == 1) xsum[[1]]$model <- object$model
		    }
                }
        class(xsum) <- "summary.rqs"
	if(class(object)[1] == "dynrqs") 
	    class(xsum) <- c("summary.dynrqs", "summary.rqs")
        xsum
        }
"logLik.rq" <- function(object,  ...){
        n <- length(object$residuals)
        p <- length(object$coefficients)
	pen <- (length(object$lambda) > 0)
	tau <- object$tau
        fid <- object$rho
        val <- n * (log(tau * (1-tau)) - 1 - log(fid/n))
        attr(val,"n") <- n
	if(pen){
	   if(!hasArg(edfThresh)) edfThresh <- 0.0001
           attr(val,"df") <- sum(abs(object$coefficients) > edfThresh)
	  }
	else  attr(val,"df") <- p
        class(val) <- "logLik"
        val
        }
"logLik.rqs" <- function(object, ...){
        n <- nrow(object$residuals)
        p <- nrow(object$coefficients)
	pen <- (length(object$lambda) > 0)
	tau <- object$tau
        fid <- object$rho
        val <- n * (log(tau * (1-tau)) - 1 - log(fid/n))
        attr(val,"n") <- n
	if(pen){
	   if(!hasArg(edfThresh)) edfThresh <- 0.0001
           attr(val,"df") <- sum(abs(object$coefficients) > edfThresh)
	  }
	else  attr(val,"df") <- p
        class(val) <- "logLik"
        val
        }
"AIC.rq" <- function(object, ... , k = 2){
        v <- logLik(object)
        if(k <= 0)
                k <- log(attr(v,"n"))
        val <- AIC(v, k = k)
        attr(val,"edf") <- attr(v,"df")
        val
        }
"extractAIC.rq"  <- function(fit, scale, k=2, ...){
aic <- AIC(fit,k)
edf <- attr(aic, "edf")
c(edf, aic)
}

"AIC.rqs" <- function(object, ... , k = 2){
        v <- logLik(object)
        if(k <= 0)
                k <- log(attr(v,"n"))
        val <- AIC(v, k = k)
        attr(val,"edf") <- attr(v,"df")
        val
        }


"summary.rq" <-
# This function provides  methods for summarizing the output of the
# rq command. In this instance, "summarizing" means essentially provision
# of either standard errors, or confidence intervals for the rq coefficents.
# Since the preferred method for confidence intervals is currently the
# rank inversion method available directly from rq() by setting ci=TRUE, with br=TRUE.
# these summary methods are intended primarily for comparison purposes
# and for use on large problems where the parametric programming methods
# of rank inversion are prohibitively memory/time consuming.  Eventually
# iterative versions of rank inversion should be developed that would
# employ the Frisch-Newton approach.
#
# Object is the result of a call to rq(), and the function returns a
# table of coefficients, standard errors, "t-statistics", and p-values, and, if
# covariance=TRUE a structure describing the covariance matrix of the coefficients,
# i.e. the components of the Huber sandwich.
#
# There are five options for "se":
#
#	1.  "rank" strictly speaking this doesn't produce a "standard error"
#		at all instead it produces a coefficient table with confidence
#		intervals for the coefficients based on inversion of the
#		rank test described in GJKP and Koenker (1994).
#	2.  "iid" which presumes that the errors are iid and computes
#		an estimate of the asymptotic covariance matrix as in KB(1978).
#	3.  "nid" which presumes local (in tau) linearity (in x) of the
#		the conditional quantile functions and computes a Huber
#		sandwich estimate using a local estimate of the sparsity.
#	4.  "ker" which uses a kernel estimate of the sandwich as proposed
#		by Powell.
#	5.  "boot" which uses a bootstrap method:
#		"xy"	uses xy-pair method
#		"wxy"	uses weighted (generalized) method
#		"pwy"	uses the parzen-wei-ying method
#		"mcmb"	uses the Markov chain marginal bootstrap method
#
#
function (object, se = NULL, covariance = FALSE, hs = TRUE, ...)
{
    if(object$method == "lasso")
         stop("no inference for lasso'd rq fitting: try rqss (if brave)")
    mt <- terms(object)
    m <- model.frame(object)
    y <- model.response(m)
    x <- model.matrix(mt, m, contrasts = object$contrasts)
    wt <- as.vector(model.weights(object$model))
    tau <- object$tau
    eps <- .Machine$double.eps^(1/2)
    coef <- coefficients(object)
    if (is.matrix(coef))
        coef <- coef[, 1]
    vnames <- dimnames(x)[[2]]
    resid <- object$residuals
    n <- length(resid)
    p <- length(coef)
    rdf <- n - p
    if (!is.null(wt)) {
        resid <- resid * wt
        x <- x * wt
        y <- y * wt
    }
    if (is.null(se)) {
        if (n < 1001 & covariance == FALSE)
            se <- "rank"
        else se <- "nid"
    }
    if (se == "rank") {
        f <- rq.fit.br(x, y, tau = tau, ci = TRUE, ...)
    }
    if (se == "iid") {
        xxinv <- diag(p)
        xxinv <- backsolve(qr(x)$qr[1:p, 1:p,drop=FALSE], xxinv)
        xxinv <- xxinv %*% t(xxinv)
        pz <- sum(abs(resid) < eps)
        h <- max(p + 1, ceiling(n * bandwidth.rq(tau, n, hs = hs)))
        ir <- (pz + 1):(h + pz + 1)
        ord.resid <- sort(resid[order(abs(resid))][ir])
        xt <- ir/(n - p)
        sparsity <- rq(ord.resid ~ xt)$coef[2]
        cov <- sparsity^2 * xxinv * tau * (1 - tau)
        scale <- 1/sparsity
        serr <- sqrt(diag(cov))
    }
    else if (se == "nid") {
        h <- bandwidth.rq(tau, n, hs = hs)
        if (tau + h > 1)
            stop("tau + h > 1:  error in summary.rq")
        if (tau - h < 0)
            stop("tau - h < 0:  error in summary.rq")
        bhi <- rq.fit.fnb(x, y, tau = tau + h)$coef
        blo <- rq.fit.fnb(x, y, tau = tau - h)$coef
        dyhat <- x %*% (bhi - blo)
        if (any(dyhat <= 0))
            warning(paste(sum(dyhat <= 0), "non-positive fis"))
        f <- pmax(0, (2 * h)/(dyhat - eps))
        fxxinv <- diag(p)
        fxxinv <- backsolve(qr(sqrt(f) * x)$qr[1:p, 1:p,drop=FALSE], fxxinv)
        fxxinv <- fxxinv %*% t(fxxinv)
        cov <- tau * (1 - tau) * fxxinv %*% crossprod(x) %*%
            fxxinv
        scale <- mean(f)
        serr <- sqrt(diag(cov))
    }
    else if (se == "ker") {
        h <- bandwidth.rq(tau, n, hs = hs)
        if (tau + h > 1)
            stop("tau + h > 1:  error in summary.rq")
        if (tau - h < 0)
            stop("tau - h < 0:  error in summary.rq")
        uhat <- c(y - x %*% coef)
        h <- (qnorm(tau + h) - qnorm(tau - h))*
		min(sqrt(var(uhat)), ( quantile(uhat,.75)- quantile(uhat, .25))/1.34 )
        f <- dnorm(uhat/h)/h
        fxxinv <- diag(p)
        fxxinv <- backsolve(qr(sqrt(f) * x)$qr[1:p, 1:p,drop=FALSE], fxxinv)
        fxxinv <- fxxinv %*% t(fxxinv)
        cov <- tau * (1 - tau) * fxxinv %*% crossprod(x) %*%
            fxxinv
        scale <- mean(f)
        serr <- sqrt(diag(cov))
    }
    else if (se == "boot") {
        B <- boot.rq(x, y, tau, ...)
        cov <- cov(B)
        serr <- sqrt(diag(cov))
    }
    if( se == "rank"){
	coef <- f$coef
	}
    else {
    	coef <- array(coef, c(p, 4))
    	dimnames(coef) <- list(vnames, c("Value", "Std. Error", "t value",
             "Pr(>|t|)"))
    	coef[, 2] <- serr
    	coef[, 3] <- coef[, 1]/coef[, 2]
    	coef[, 4] <- if (rdf > 0)
			2 * (1 - pt(abs(coef[, 3]), rdf))
    		     else NA
	}
    object <- object[c("call", "terms")]
    if (covariance == TRUE) {
        object$cov <- cov
        if(se == "iid") object$scale <- scale
        if(se %in% c("nid", "ker")) {
            object$Hinv <- fxxinv
            object$J <- crossprod(x)
            object$scale <- scale
        }
        else if (se == "boot") {
            object$B <- B
        }
    }
    object$coefficients <- coef
    object$residuals <- resid
    object$rdf <- rdf
    object$tau <- tau
    class(object) <- "summary.rq"
    object
}

akj <- function(x,
                z = seq(min(x), max(x), length = 2 * length(x)),
                p = rep(1/length(x), length(x)),
                h = -1, alpha = 0.5, kappa = 0.9, iker1 = 0)
{
    nx <- length(x)
    stopifnot(is.numeric(x),
              length(p) == nx,
	      any((iker1 <- as.integer(iker1)) == 0:1))
    nz <- length(z)
    if(is.unsorted(x))
	x <- sort(x)
    .Fortran("akj",
	     as.double(x),
	     as.double(z),
	     as.double(p),
	     iker1,
	     dens  = double(nz),
	     psi   = double(nz),
	     score = double(nz),
	     as.integer(nx),
	     as.integer(nz),
	     h = as.double(h),
	     as.double(alpha),
	     as.double(kappa),
	     double(nx),
	     PACKAGE = "quantreg")[c("dens", "psi", "score", "h")]
}

"lm.fit.recursive" <-
function(X, y, int = TRUE)
{
	if(int)
		X <- cbind(1, X)
	p <- ncol(X)
	n <- nrow(X)
	D <- qr(X[1:p,  ])
	A <- qr.coef(D, diag(p))
	A[is.na(A)] <- 0
	A <- crossprod(t(A))
	Ax <- rep(0, p)
	b <- matrix(0, p, n)
	b[, p] <- qr.coef(D, y[1:p])
	b[is.na(b)] <- 0
	z <- .Fortran( "rls",
		as.integer(n),
		as.integer(p),
		as.double(t(X)),
		as.double(y),
		b = as.double(b),
		as.double(A),
		as.double(Ax),
		PACKAGE = "quantreg")
	bhat <- matrix(z$b, p, n)
	return(bhat)
}

"rq.fit.hogg" <-
function (x, y, taus = c(.1,.3,.5), weights = c(.7,.2,.1),  
	R= NULL, r = NULL, beta = 0.99995, eps = 1e-06) 
{
    n <- length(y)
    n2 <- nrow(R)
    m <- length(taus)
    p <- ncol(x)+m
    if (n != nrow(x)) 
        stop("x and y don't match n")
    if (m != length(weights)) 
        stop("taus and weights differ in length")
    if (any(taus < eps) || any(taus > 1 - eps)) 
        stop("taus outside (0,1)")
    W <- diag(weights)
    if(m == 1) W <- weights
    x <- as.matrix(x)
    X <- cbind(kronecker(W,rep(1,n)),kronecker(weights,x))
    y <- kronecker(weights,y)
    rhs <- c(weights*(1 - taus)*n, sum(weights*(1-taus)) * apply(x, 2, sum))
    if(n2!=length(r))
	stop("R and r of incompatible dimension")
    if(ncol(R)!=p)
	stop("R and X of incompatible dimension")
    d <- rep(1, m*n)
    u <- rep(1, m*n)
    if(length(r)){
       wn1 <- rep(0, 10 * m*n)
       wn1[1:(m*n)] <- .5
       wn2 <- rep(0,6*n2)
       wn2[1:n2] <- 1 
       z <- .Fortran("rqfnc", as.integer(m*n), as.integer(n2), as.integer(p), 
           a1 = as.double(t(as.matrix(X))), c1 = as.double(-y), 
           a2 = as.double(t(as.matrix(R))), c2 = as.double(-r), 
           rhs = as.double(rhs), d1 = double(m*n), d2 = double(n2), 
           as.double(u), beta = as.double(beta), eps = as.double(eps), 
           wn1 = as.double(wn1), wn2 = as.double(wn2), wp = double((p + 3) * p), 
	   it.count = integer(3), info = integer(1), PACKAGE = "quantreg")
	}
    else{
	wn <- rep(0, 10 * m*n)
    	wn[1:(m*n)] <- .5
    	z <- .Fortran("rqfnb", as.integer(m*n), as.integer(p), a = as.double(t(as.matrix(X))), 
		c = as.double(-y), rhs = as.double(rhs), d = as.double(d), as.double(u), 
		beta = as.double(beta), eps = as.double(eps), wn = as.double(wn), 
		wp = double((p + 3) * p), it.count = integer(2), info = integer(1),
		PACKAGE = "quantreg")
	}
    if (z$info != 0) 
        warning(paste("Info = ", z$info, "in stepy: singular design: iterations ", z$it.count[1]))
    coefficients <- -z$wp[1:p]
    if(any(is.na(coefficients)))stop("NA coefs:  infeasible problem?")
    list(coefficients = coefficients, nit = z$it.count, flag = z$info)
}
