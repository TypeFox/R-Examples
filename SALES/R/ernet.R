ernet <- function(x, y, nlambda = 100L, method = "er", 
                  lambda.factor = ifelse(nobs < nvars, 1e-02, 1e-04), 
                  lambda = NULL, lambda2 = 0, pf = rep(1, nvars), 
                  pf2 = rep(1, nvars), exclude, dfmax = nvars + 1, 
                  pmax = min(dfmax * 1.2, nvars), standardize = TRUE, 
                  intercept = TRUE, eps = 1e-08, maxit = 1000000L, 
                  tau = 0.5) {
    #################################################################################
    ## data setup
    method <- match.arg(method)
    this.call <- match.call()
    y <- drop(y)
    x <- as.matrix(x)
    np <- dim(x)
    nobs <- as.integer(np[1])
    nvars <- as.integer(np[2])
    vnames <- colnames(x)
    if (is.null(vnames)) vnames <- paste("V", seq(nvars), sep = "")
    if (NROW(y) != nobs) stop("x and y have different number of observations")
    if (NCOL(y) > 1L) stop("Multivariate response is not supported now")
    #################################################################################
    ## parameter setup
    if (length(pf) != nvars) 
      stop("Size of L1 penalty factors does not match the number of input variables")
    if (length(pf2) != nvars) 
      stop("Size of L2 penalty factors does not match the number of input variables")
    if (lambda2 < 0) stop("lambda2 should be non-negative")
    maxit <- as.integer(maxit)
    lam2 <- as.double(lambda2)
    pf <- as.double(pf)
    pf2 <- as.double(pf2)
    isd <- as.integer(standardize)
    intr <- as.integer(intercept)
    eps <- as.double(eps)
    dfmax <- as.integer(dfmax)
    pmax <- as.integer(pmax)
    if (!missing(exclude)) {
      jd <- match(exclude, seq(nvars), 0)
      if (!all(jd > 0)) stop("Some excluded variables out of range")
      jd <- as.integer(c(length(jd), jd))
    } else jd <- as.integer(0)
    #################################################################################
    ## lambda setup
    nlam <- as.integer(nlambda)
    if (is.null(lambda)) {
      if (lambda.factor >= 1) stop("lambda.factor should be less than 1")
      flmin <- as.double(lambda.factor)
      ulam <- double(1)
    } else {
        flmin <- as.double(1) # flmin = 1 if user defines lambda
        if (any(lambda < 0)) stop("lambdas should be non-negative")
        ulam <- as.double(rev(sort(lambda)))
        nlam <- as.integer(length(lambda))
    }
    #################################################################################
    fit <- alspath(x, y, nlam, flmin, ulam, isd, intr, eps, dfmax, pmax, jd, 
                pf, pf2, maxit, lam2, tau, nobs, nvars, vnames)
    if (is.null(lambda)) 
        fit$lambda <- lamfix(fit$lambda)
    fit$call <- this.call
    #################################################################################
    class(fit) <- c("ernet", class(fit))
    fit
} 
