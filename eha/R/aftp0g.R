aftp0g <- function(printlevel, ns, nn, id,
                   strata, Y, X, offset, means, param){
    pfixed <- FALSE # Never fix 'shape' with Gompertz (why?)
    ## Gompertz aft-reg!
    Fmin <- function(beta){

        fit <- .C("aftregGomp",
                  as.integer(printlevel),
                  ##
                  as.integer(ns), # No. of strata
                  as.integer(nn),
                  as.integer(ncov),
                  as.integer(bdim),
                  ##
                  as.integer(id),
                  as.integer(strata - 1), # 0, ..., (ns-1); C-style!
                  as.double(Y[, 1]),  ## 'enter'
                  as.double(Y[, 2]),  ## 'exit'
                  as.integer(Y[, 3]), ## 'event'
                  ##
                  ##as.double(t(scale(X, center = TRUE, scale = FALSE))),
                  as.double(t(X)), # NOTE; scaling already done!
                  ## NOTE transpose!
                  as.double(offset),
                  as.integer(dis),     # baseline distribution
                  as.double(beta),
                                        # results -->
                  loglik = double(1),  # function value at beta
                  fail = integer(1), # = 0: No failure
                  ## DUP = FALSE,
                  PACKAGE = "eha"
                  )
        if (fit$fail) stop("Error in likelihood calculation")
        return(fit$loglik)
    }
    ncov <- NCOL(X)
    dis <- 4 ## Just a joke to make Gomp happy
    ## Start values (no covariates):
    bdim <- 2 * ns
    ncov.save <- ncov
    ncov <- 0
    beta <- numeric(bdim)
    for (i in seq_len(ns)){
        enter <- Y[strata == i, 1]
        exit <-  Y[strata == i, 2]
        event <-  Y[strata == i, 3]
        
        #alpha <- log(max(exit))
        beta[(2 * i - 1):(2 * i)] <- gompstartAft(enter, exit, event)
        ##beta[2 * i - 1] <- log(sum(Y[, 2] - Y[, 1]) / sum(Y[, 3]))
        ##beta[2*i] <- 0
    }
    ## This is not strictly becessary, but:
    res0 <- optim(beta, Fmin, method = "BFGS", hessian = FALSE)
    loglik.start <- -res0$value
    ##cat("\nloglik.start = ", loglik.start, "\n")
    ##cat("beta = ", res0$par, "\n\n")
    ##cat("convergence = ", res0$convergence, "\n\n")
    ## And now the real thing:
    ncov <- ncov.save
    bdim <- ncov + 2 * ns
    beta <- c(rep(0, ncov), beta)
    res <- optim(beta, Fmin, method = "BFGS",
                 control = list(trace = as.integer(printlevel)),
                 hessian = TRUE)
    ##cat("\nloglik = ", -res$value, "\n")
    fit <- list(beta = res$par, loglik = c(loglik.start, -res$value))

    fit$var <- try(solve(res$hessian))
    fit$fail <- res$convergence != 0 # Optimistic?
        
    fit$coef.names <- names(fit$beta)

    fit$ncov <- ncov
    fit
}
