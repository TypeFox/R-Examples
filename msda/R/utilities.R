# KKT check the output theta matrix
kktchk <- function(obj, pf, thr) {
    lambda <- obj$lambda
    nk <- nrow(obj$delta)
    nvars <- ncol(obj$sigma)
    for (l in 1:length(lambda)) {
        cat("now checking lambda ", l, "\n")
        theta <- t(obj$theta[[l]])
        thetaInner <- apply(theta, 2, crossprod)
        thetaNorm <- sqrt(thetaInner)
        
        sdiag <- diag(obj$sigma)
        sjj <- t(replicate(nk, sdiag))
        thetaTilda <- (obj$delta - theta %*% obj$sigma + sjj * theta)/sjj
        thetaDif <- theta - thetaTilda
        
        for (j in 1:nvars) {
            los <- lambda[l] * pf[j]/obj$sigma[j, j]
            if (thetaNorm[j] == 0) {
                dif_norm <- sqrt(crossprod(thetaDif[, j]))
                tmp <- dif_norm - los
                if (tmp > 0) 
                  cat("violated at t > 0", tmp, "\n")
            } else {
                tmp3 <- thetaDif[, j] + los * theta[, j]/thetaNorm[j]
                if (any(abs(tmp3) > thr)) 
                  cat("violated at t = 0", tmp3, "\n")
            }
        }
    }
}

# extract fortran outputs and format it into sparse matries
formatoutput <- function(fit, maxit, pmax, nvars, vnames, nk) {
    nalam <- fit$nalam
    ntheta <- fit$ntheta[seq(nalam)]
    nthetamax <- max(ntheta)
    lam <- fit$alam[seq(nalam)]
	obj <- fit$obj[seq(nalam)]
    stepnames <- paste("s", seq(nalam) - 1, sep = "")
    resnames <- paste("delta", seq(nk), sep = "")
    
    errmsg <- err(fit$jerr, maxit, pmax)  ### error messages from fortran
    switch(paste(errmsg$n), `1` = stop(errmsg$msg, call. = FALSE), `-1` = cat(errmsg$msg))
    
    dd <- c(nvars, nk)
    df <- rep(0, nalam)
    if (nthetamax > 0) {
        ja <- fit$itheta[seq(nthetamax)]
        oja <- order(ja)
        ja <- rep(ja[oja], nk)
        itheta <- cumsum(c(1, rep(nthetamax, nk)))
        pos <- rep(1:nalam, each = nk * pmax)
        theta <- split(fit$theta[seq(nk * pmax * nalam)], pos)
        for (l in 1:nalam) {
            theta[[l]] <- matrix(theta[[l]], pmax, nk, byrow = TRUE)[seq(nthetamax), 
                , drop = FALSE]
            df[l] <- sum(rowSums(abs(theta[[l]])) != 0)
            theta[[l]] <- new("dgCMatrix", Dim = dd, Dimnames = list(vnames, 
                resnames), x = as.vector(theta[[l]][oja, ]), p = as.integer(itheta - 
                1), i = as.integer(ja - 1))
        }
    } else {
        theta <- list()
        for (l in 1:nalam) {
            theta[[l]] <- zeromat(nvars, nk, vnames, resnames)
        }
        df <- rep(0, nalam)
    }
    list(theta = theta, df = df, dim = dd, lambda = lam, obj = obj)
}


# generate sigma, delta and mu from x, y
msda.prep <- function(x, y) {
    # data setup
    x <- as.matrix(x)
    y <- drop(y)
    nclass <- as.integer(length(unique(y)))
    prior <- rep(0, nclass)
    for (k in 1:nclass) {
        prior[k] <- mean(y == k)
    }
    nvars <- as.integer(ncol(x))
    nobs <- nrow(x)
    nres <- length(y)
    if (nres != nobs) 
        stop("x and y have different number of observations")
    # prepare sigma and delta
    mu <- matrix(0, nvars, nclass)
    sigma <- matrix(0, nvars, nvars)
    for (i in 1:nclass) {
        mu[, i] <- apply(x[y == i, ], 2, mean)
        sigma <- sigma + (sum(y == i) - 1) * cov(x[y == i, ])
    }
    sigma <- sigma/(nobs - nclass)
    delta <- sweep(mu[, -1], 1, mu[, 1], "-")
    delta <- t(delta)
    outlist <- list(sigma = sigma, delta = delta, mu = mu, prior = prior)
    outlist
}

err <- function(n, maxit, pmax) {
    if (n == 0) 
        msg <- ""
    if (n > 0) {
        # fatal error
        if (n < 7777) 
            msg <- "Memory allocation error; contact package maintainer"
        if (n == 10000) 
            msg <- "All penalty factors are <= 0"
        n <- 1
        msg <- paste("in the fortran code -", msg)
    }
    if (n < 0) {
        # non fatal error
        if (n > -10000) 
            msg <- paste("Convergence for ", -n, "th lambda value not reached after maxit=", 
                maxit, " iterations; solutions for larger lambdas returned.\n", 
                sep = "")
        if (n < -10000) 
            msg <- paste("Number of nonzero coefficients along the path exceeds pmax=", 
                pmax, " at ", -n - 10000, "th lambda value; solutions for larger lambdas returned.\n", 
                sep = "")
        if (n < -20000) 
            msg <- paste("Number of nonzero coefficients along the path exceeds dfmax=", 
                pmax, " at ", -n - 20000, "th lambda value; solutions for larger lambdas returned.\n", 
                sep = "")
        n <- -1
    }
    list(n = n, msg = msg)
}

zeromat <- function(nvars, nalam, vnames, stepnames) {
    ca <- rep(0, nalam)
    ia <- seq(nalam + 1)
    ja <- rep(1, nalam)
    dd <- c(nvars, nalam)
    new("dgCMatrix", Dim = dd, Dimnames = list(vnames, stepnames), x = as.vector(ca), 
        p = as.integer(ia - 1), i = as.integer(ja - 1))
}

lamfix <- function(lam) {
    llam <- log(lam)
    lam[1] <- exp(2 * llam[2] - llam[3])
    lam
}
