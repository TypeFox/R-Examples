nnnegbin <- function(y, x, standard, offset, start, control = list())
{
	control <- do.call("addreg.control", control)
	x <- as.matrix(x)
	xnames <- dimnames(x)[[2L]]
	ynames <- if (is.matrix(y))
		rownames(y)
	else names(y)
	
	if (any(x < 0)) stop("x must be non-negative")
	if (any(apply(x, 2, function(col) all(col == 0)))) stop("x contains column with all 0")
	
	nvars <- ncol(x)
	nobs <- NROW(y)
	
	fam <- negbin1(link = identity, phi = NA)
	eval(fam$initialize)
	
	weights <- rep(1, nobs)
	if (is.null(standard)) standard <- rep.int(1, nobs)
	if (!is.null(offset))
		stop("offset is not currently supported for negbin1 family")
		
	converged <- FALSE
	
	if (!is.null(start)) {
		if (length(start) != nvars + 1)
			stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s",
				nvars + 1, paste(c(deparse(xnames), "scale"), collapse = ", ")),
				domain = NA)
		else if (any(start[-length(start)] <= control$bound.tol))
			stop("'start' is on our outside the boundary of the parameter space (consider 'bound.tol')",
				domain = NA)
		else if (start[length(start)] <= 1)
			stop("'start' value for scale must be > 1 for the negbin1 family", domain = NA)
		else {
			coefold.mu <- start[-length(start)]
			coefold.phi <- start[length(start)] - 1
		}
	} else {
		simple <- mean(y / standard) / colMeans(x) + 2*control$bound.tol
		trymat <- tryCatch(as.vector(solve(t(x)%*%x) %*% t(x) %*% (y)) + 2*control$bound.tol, error = function(e) NULL)
		if (is.null(trymat) || any(trymat < control$bound.tol)) coefold.mu <- simple
		else coefold.mu <- trymat
		coefold.phi <- 0.5
	}
	
	coefold.r <- coefold.mu / coefold.phi
	coefold.p <- coefold.phi / (1 + coefold.phi)
	r.fits <- standard * drop(x %*% coefold.r)
	
	fam <- negbin1(link = identity, phi = coefold.phi)
	mu.eta <- fam$mu.eta
	linkinv <- fam$linkinv
	dev.resids <- fam$dev.resids
	aic <- fam$aic
	
	eta <- drop(x %*% coefold.mu)
	mu <- standard * linkinv(eta)
	dev.old <- sum(dev.resids(y, mu, weights))
	
	ll.old <- -aic(y, nobs, mu, weights, dev.old)/2
	
	if (control$trace) cat("Log-likelihood =", ll.old, "Iterations - 0\n")
	
	dbetabinom <- function(x, size, alpha, beta) {
		exp(lchoose(size, x) + lbeta(alpha + x, size + beta - x) - lbeta(alpha, beta))
	}
	
	digammadiff <- function(y, r) {
		k <- seq_len(y) - 1
		sum(1 / (r + k))
	}
    
    digammaexp <- function(y, r, r.part, r.full) {
        k <- seq_len(y)
        rdiff <- cumsum(1 / (r + k - 1))
        prob <- sapply(k, dbetabinom, size = y, alpha = r.part, beta = r.full - r.part)
        sum(rdiff * prob)
    }
	
	estep <- function(theta, x, y, std, theta.old, r.fits.old) {
		res <- rep(0, length(y))
		r.part.old <- std * theta.old * x
		chk1 <- (y == 0)
		chk2 <- (r.part.old == r.fits.old)
		if(sum(!chk1 & chk2) > 0)
			res[!chk1 & chk2] <- mapply(digammadiff, y = y[!chk1 & chk2],
										r = (std * theta * x)[!chk1 & chk2])
		if (sum(!chk1 & !chk2) > 0)
			res[!chk1 & !chk2] <- mapply(digammaexp, y = y[!chk1 & !chk2],
                                        r = (std * theta * x)[!chk1 & !chk2],
                                        r.part = r.part.old[!chk1 & !chk2],
                                        r.full = r.fits.old[!chk1 & !chk2])
		res
	}
		
	score <- function(theta, x, y, std, theta.old, r.fits.old, p) {
		sum(std[x > 0] * x[x > 0] * (estep(theta, x[x > 0], y[x > 0], std[x > 0], theta.old,
			r.fits.old[x > 0]) + log(1 - p)))
	}
	
	for (iter in 1L:control$maxit) {
		coefmax <- drop(coefold.r * (t(ifelse(r.fits == 0, 0, standard * y / r.fits)) %*% x)
					/ (log(1 / (1 - coefold.p)) * t(standard) %*% x))
		coefnew.r <- coefold.r
		for (j in 1L:nvars) {
			if (coefold.r[j] > control$bound.tol) {
				score.0 <- score(control$bound.tol / 2, x = x[,j], y = y, std = standard, 
					theta.old = coefold.r[j], r.fits.old = r.fits, p = coefold.p)
				score.max <- score(coefmax[j], x = x[,j], y = y, std = standard,
					theta.old = coefold.r[j], r.fits.old = r.fits, p = coefold.p)
				if (score.0 <= control$epsilon)
					coefnew.r[j] <- 0
				else if (score.max >= 0)
					coefnew.r[j] <- coefmax[j]
				else
					coefnew.r[j] <- uniroot(score, interval = c(control$bound.tol / 2, coefmax[j]),
						x = x[,j], y = y, std = standard, theta.old = coefold.r[j],
						r.fits.old = r.fits, p = coefold.p, f.lower = score.0, f.upper = score.max,
						tol = control$epsilon * 1e-2)$root
			}
		}
		
		r.fits <- standard * drop(x %*% coefnew.r)
		coefnew.p <- sum(y) / (sum(y) + sum(r.fits))
		
		coefnew.phi <- coefnew.p / (1 - coefnew.p)
		coefnew.mu <- coefnew.phi * coefnew.r
		
		fam <- negbin1(link = identity, phi = coefnew.phi)
		mu.eta <- fam$mu.eta
		linkinv <- fam$linkinv
		dev.resids <- fam$dev.resids
		aic <- fam$aic
		
		eta <- drop(x %*% coefnew.mu)
		mu <- standard * linkinv(eta)
		dev.new <- sum(dev.resids(y, mu, weights))
		ll.new <- -aic(y, nobs, mu, weights, dev.new)/2
		
		if (control$trace) cat("Log-likelihood =", ll.new, "Iterations -", iter, "\n")
		
        if (conv.test(coefold.mu, coefnew.mu, control$epsilon)
            && conv.test(coefold.phi, coefnew.phi, control$epsilon)) {
			converged = TRUE
			break
		}
		
		coefold.r <- coefnew.r
		coefold.p <- coefnew.p
		coefold.mu <- coefnew.mu
		coefold.phi <- coefnew.phi
	}
	
	residuals <- (y - mu) / mu.eta(eta)
	
	names(coefnew.mu) <- xnames
	names(residuals) <- names(mu) <- names(eta) <- names(y) <- ynames
	
	aic.model <- aic(y, nobs, mu, weights, dev.new) + 2 * (nvars + 1)
	aic.c.model <- aic.model + 2 * (nvars + 1) * (nvars + 2) / (nobs - nvars - 1)
	
	wtdmu <- standard * rep(sum(weights * y / standard) / sum(weights), nobs)
	nulldev <- sum(dev.resids(y, wtdmu, weights))
	nulldf <- nobs - 1
	resdf <- nobs - nvars - 1
	
	boundary <- any(coefnew.r < control$bound.tol)
	
	list(coefficients = coefnew.mu, scale = 1 + coefnew.phi, residuals = residuals, fitted.values = mu,
		rank = nvars + 1, family = fam, linear.predictors = eta, deviance = dev.new, aic = aic.model,
		aic.c = aic.c.model, null.deviance = nulldev, iter = iter, weights = weights, prior.weights = weights,
		standard = standard, df.residual = resdf, df.null = nulldf, y = y, converged = converged, boundary = boundary,
		loglik = ll.new, nn.design = x)
}