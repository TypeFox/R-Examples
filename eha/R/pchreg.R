pchreg <- function(X, Y, cuts, offset, init, control, center){
    ## Piecewise constant (pch) proportional hazards:
    ## h(t; beta) = 1/b_i * exp(t / b_i) * exp(x*beta),
    ## i = 1, ..., (No. of constant intervals).
    ## NOTE: Y is nn x 3!
    
    if (!is.matrix(X)) X <- matrix(X, ncol = 1)
    nn <- NROW(X)
    ncov <- NCOL(X)
    if (NROW(Y) != nn) stop("Y NROW error")
    if (NCOL(Y) != 3) stop("Y NCOL error", )
    if (ncov){
        wts <- Y[, 2] - Y[, 1]
        w.means <- apply(X, 2, weighted.mean, w = wts)
        means <- apply(X, 2, mean)
        ## Changed in 2.2-1 (23 March, 2013); no centering:
        ##if (center){
        ##    for (i in 1:ncov){
        ##        X[, i] <- X[, i] - means[i]
        ##    }
        ##}
        ## Will be dealt with (internal centering) later.
    }
    strata <- rep(1, nn)
    ns <- 1

    if (missing(offset) || is.null(offset))
        offset <- rep(0, nn)

    if (missing(init) || is.null(init))
        init <- rep(0, ncov)
    if (length(init) != ncov) stop("Error in init")

    printlevel <- control$trace
    iter <- control$maxiter


    bdim <- ncov ## + 2

    ## First, fit the 'null' model: NO!

    ##ncov.save <- ncov
    ##ncov <- 0

    split <- SurvSplit(Y, cuts)
    event <- split$Y[, 3]
    Y <- split$Y
    X <- X[split$idx, ,drop = FALSE]
    offset <- offset[split$idx] + log(Y[, 2] - Y[, 1])
    res0 <- glmmboot(event ~ offset(offset), cluster = split$ivl,
                     family = "poisson")

    ## Done; now the real thing: (The above is really unnecessary? No!!)
    ##ncov <- ncov.save
    
    ##if (ncov)
      ##  beta[1:ncov] <- init  # Start values
    res <- glmmbootFit(X, event, cluster = split$ivl,
                       offset = offset, family = poisson())
    coefficients <- res$coefficients
    coef.names <- colnames(X)
    names(coefficients) <- coef.names
    fit <- list(coefficients = coefficients,
                loglik = c(res0$logLik, res$logLik),
                hazards = exp(res$frail)
                )
    fit$pfixed <- TRUE
    fit$var <- res$variance
        ##tryCatch(solve(-res$hessian), error = function(e) e)
    if (is.matrix(fit$var)){
        colnames(fit$var) <- coef.names
        rownames(fit$var) <- coef.names
    }
    ##fit$hessian <- res$hessian
    ##if (is.matrix(fit$hessian)){
    ##    colnames(fit$hessian) <- coef.names
    ##    rownames(fit$hessian) <- coef.names
    ##}
    
    ##fit$n.strata <- ns
    fit$df <- ncov
    fit$fail <- FALSE # Optimist!
    fit$cuts <- cuts
    fit
}
                      
