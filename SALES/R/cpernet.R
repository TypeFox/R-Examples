cpernet <- function(x, y, w = 1.0, nlambda = 100L, method = "cper", 
                    lambda.factor = ifelse(2 * nobs < nvars, 1e-02, 1e-04), 
                    lambda = NULL, lambda2 = 0, pf.mean = rep(1, nvars), 
                    pf2.mean = rep(1, nvars), pf.scale = rep(1, nvars),
                    pf2.scale = rep(1, nvars), exclude, dfmax = nvars + 1, 
                    pmax = min(dfmax * 1.2, nvars), standardize = TRUE, 
                    intercept = TRUE, eps = 1e-08, maxit = 1000000L, 
                    tau = 0.80) {
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
    if (length(pf.mean) != nvars) 
      stop("Size of L1 penalty factors for the mean does not match the number of input variables")
	  if (length(pf2.mean) != nvars) 
      stop("Size of L2 penalty factors for the mean does not match the number of input variables")
    if (length(pf.scale) != nvars) 
      stop("Size of L1 penalty factors for the scale does not match the number of input variables")
    if (length(pf2.scale) != nvars) 
      stop("Size of L2 penalty factors for the scale does not match the number of input variables")
    if (lambda2 < 0) {
       warning("lambda2 < 0; set to zero...")
       lambda2 <- 0
    }
    maxit <- as.integer(maxit)
    lam2 <- as.double(lambda2)
    pfmean <- as.double(pf.mean)
    pf2mean <- as.double(pf2.mean)
    pfscale <- as.double(pf.scale)
    pf2scale <- as.double(pf2.scale)
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
    fit <- cpalspath(x, y, w, nlam, flmin, ulam, isd, intr, eps, dfmax, pmax, jd, 
                     pfmean, pf2mean, pfscale, pf2scale, maxit, lam2, tau, 
                     nobs, nvars, vnames)
    if (is.null(lambda)) fit$lambda <- lamfix(fit$lambda)
    fit$call <- this.call
    #################################################################################
    class(fit) <- c("cpernet", class(fit))
    fit
} 
