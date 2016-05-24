nplbin <- function(y, x, offset, start, control = list())
{
    control <- do.call("logbin.control", control)
    x <- as.matrix(x)
    xnames <- dimnames(x)[[2L]]
    ynames <- if (is.matrix(y))
        rownames(y)
    else names(y)
    
    if (any(x < 0)) stop("x must be non-negative")
    if (any(apply(x, 2, function(col) all(col==0)))) stop("x contains column with all 0")
    
    nvars <- ncol(x)
    nobs <- NROW(y)
    
    n <- weights <- rep(1, nobs)
    if (is.null(offset)) offset <- rep.int(0, nobs)
    if (any(offset > 0))
        stop("offset must be non-positive")
    
    fam <- binomial(link = log)
    eval(fam$initialize)
    
    mu.eta <- fam$mu.eta
    linkinv <- fam$linkinv
    dev.resids <- fam$dev.resids
    aic <- fam$aic
    loglik <- function(y, mu, n) sum(lchoose(n,y) + y*log(mu) + (n-y)*log(1-mu))
    
    y1 <- round(n*y)
    y2 <- round(n*(1-y))
    
    x.mins <- apply(x, 2, function(t) min(t[t > 0]))
    x.scale <- 1/x.mins
    x1 <- sweep(x, 2, FUN = "*", x.scale)
    
    converged <- FALSE
    
    coefold <- if (!is.null(start)) {
        if (length(start) != nvars)
            stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s",
                    nvars, paste(deparse(xnames), collapse = ", ")),
                domain = NA)
        else if (any(start >= -control$bound.tol))
            stop("'start' is on our outside the boundary of the parameter space (consider 'bound.tol')", domain = NA)
        else start
    } else {
        simple <- log(mean(y)) / colMeans(x) - 2*control$bound.tol
        logy <- log(y)
        logy[is.infinite(logy)] <- min(logy[!is.infinite(logy)]) - 2*control$bound.tol
        trymat <- tryCatch(as.vector(solve(t(x) %*% x) %*% t(x) %*% (logy)) + 2*control$bound.tol,
                    error = function(e) NULL)
        if (is.null(trymat)) simple
        else if (any(trymat >= -control$bound.tol)) simple
        else trymat
    }
    
    coefold <- coefold / x.scale
    
    eta <- drop(x1 %*% coefold) + offset
    mu <- n * linkinv(eta)
    dev.old <- sum(dev.resids(y, mu / n, n))
    
    ll.old <- loglik(y1, mu / n, n)
    
    if (control$trace) cat("Deviance =", dev.old, "Iterations - 0\n")
    
    std.div <- 1 / colSums(n * x1)
    
    for (iter in 1L:control$maxit) {
    
        estep <- y1 + y2*((matrix(linkinv(coefold), nobs, nvars, byrow = TRUE) - linkinv(eta))/(1 - linkinv(eta)))
        coefnew <- log(colSums(estep * x1) * std.div)
		coefnew[coefnew >= 0] <- -control$bound.tol / 2
        
        eta <- drop(x1 %*% coefnew) + offset
        mu <- n * linkinv(eta)
        dev.new <- sum(dev.resids(y, mu / n, n))
        
        ll.new <- loglik(y1, mu / n, n)
        
        if(control$trace) cat("Deviance =", dev.new, "Iterations -", iter, "\n")
        
        if(conv.test(coefold * x.scale, coefnew * x.scale, control$epsilon)) {
            converged = TRUE
            break
        }
        
        coefold <- coefnew
        ll.old <- ll.new
        dev.old <- dev.new
    }
    
    coefnew <- coefnew * x.scale
    
    residuals <- (y - (mu / n)) / mu.eta(eta)
    
    names(coefnew) <- xnames
    names(y) <- names(mu) <- names(eta) <- names(residuals) <- ynames
    
    aic.model <- aic(y, n, mu / n, weights, dev.new) + 2 * nvars
    aic.c.model <- aic.model + 2 * nvars * (nvars + 1) / (nobs - nvars - 1)
    
    wtdmu <- sum(n * y) / sum(n)
    nulldev <- sum(dev.resids(y, wtdmu, n))
    nulldf <- nobs - 1
    resdf <- nobs - nvars
    
    boundary <- any(coefnew > -control$bound.tol)
    
    list(coefficients = coefnew, residuals = residuals, fitted.values = mu / n, rank = nvars,
            family = fam, linear.predictors = eta, deviance = dev.new, aic = aic.model, 
            aic.c = aic.c.model, null.deviance = nulldev, iter = iter, weights = weights, 
            prior.weights = n, df.residual = resdf, df.null = nulldf, y = y, 
            converged = converged, boundary = boundary, loglik = ll.new, nn.design = x)   
}