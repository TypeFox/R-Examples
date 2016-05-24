"grouped.fit" <-
function(y, X, link, distribution, df., starts, iter){
    if(is.null(n <- nrow(X)))
        stop("`X' must be a matrix.\n")
    if(n == 0)
        stop("0 (non-NA) cases.\n")
    px <- ncol(X)
    p <- px + 1
    if(qr(X)$rank < px)
        stop("`X' is singular.\n")
    if(px == 0)
        return(list(coefficients = NULL, logLik = NULL, hessian = NULL, k = NULL, n = n))
    if(!is.matrix(y) || ncol(y) != 2)
        stop("The response must be a two-column matrix.\n")
    if(n != nrow(y))
        stop("the dimensions of `y' and `X' do not mutch.\n")
    a <- y[, 1]
    b <- y[, 2]
    if(any(a > b)) stop("The values in the response are not valid;")
    qa <- switch(link, "identity" = a, "log" = log(a), "logit" = qlogis(a))
    qb <- switch(link, "identity" = b, "log" = log(b), "logit" = qlogis(b))
    qa. <- ifelse(is.finite(qa), qa, 0.0)
    qb. <- ifelse(is.finite(qb), qb, 0.0)
    start.par <- if(missing(starts)) start.values(a, b, X) else{
        if(length(starts) != p){
            warning("The starting values and the number of parameters do not much, random starting values are used instead.\n")
            start.values(a, b, X)
        } else starts
    }
    old <- options(warn = (-1))
    on.exit(options(old))
    lL <- switch(distribution, "normal" = lL.gaussian, "t" = lL.t, "logistic" = lL.logistic)
    sc <- switch(distribution, "normal" = sc.gaussian, "t" = sc.t, "logistic" = sc.logistic)
    environment(lL) <- environment()
    environment(sc) <- environment()
    fn0 <- lL(start.par)
    gn0 <- sc(start.par)
    while(!all(is.finite(c(fn0, gn0)))){
        start.par <- c(rnorm(px), 1)
        fn0 <- lL(start.par)
        gn0 <- sc(start.par)
    }
    opt <- optim(start.par, fn = lL, gr = sc, method = "BFGS", hessian = TRUE, 
                    control = list(maxit = 250, reltol = 1e-10))
    k <- 1
    while(opt$conv != 0){
        opt <- optim(opt$par, fn = lL, gr = sc, method = "BFGS", hessian = TRUE)
        k <- k + 1
        if(k == iter) break
    }
    if(conv <- opt$conv)
        warning("The algorithm did not convergence; re-fit using different starting values.\n")
    ev <- eigen(hes <- opt$hes, TRUE, TRUE)$values
    if(!all(ev >= -1e-06 * abs(ev[1]))) 
        warning("Hessian matrix at convergence is not positive definite, unstable solution.\n")
    details <- list(X = X, y = y, convergence = conv, logLik = -opt$val, k = k, n = n, df = df., link = link, 
                        distribution = distribution, max.sc = max(abs(sc(opt$par))))
    theta <- opt$par
    fits <- switch(link, "identity" = c(X %*% theta[1:px]),
                         "log" = c(exp(X %*% theta[1:px])),
                         "logit" = c(plogis(X %*% theta[1:px])))
    dn <- colnames(X)
    xnames <- if(is.null(dn)) c("(Intercept)", paste("X", seq(1, px - 1), sep = "")) else dn
    xnames <- c(xnames, "sigma")
    names(theta) <- xnames
    dimnames(hes) <- list(xnames, xnames)
    names(fits) <- rownames(X)
    fit <- list(coefficients = theta, hessian = hes, fitted = fits, details = details)
    fit
}

