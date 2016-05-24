aftp1 <- function(printlevel, ns, nn, id,
                  strata, Y, X, offset, shape, dis, means){
    Fexpmin <- function(beta){

        fit <- .C("aftexpsup",
                  as.integer(printlevel),
                                        #
                  as.integer(ns), # No. of strata
                  as.integer(nn),
                  as.integer(ncov),
                  as.integer(bdim),
                                        #
                  as.integer(id),
                  as.integer(strata - 1), # 0, ..., (ns-1); C-style!
                  as.double(Y[, 1]),  ## 'enter'
                  as.double(Y[, 2]),  ## 'exit'
                  as.integer(Y[, 3]), ## 'event'
                                        #
                  ##as.double(t(scale(X, center = TRUE, scale = FALSE))),
                  as.double(t((X))), #NOTE: scaling already done!
                  as.double(offset),
                  as.double(shape), ## "p" (fixed)
                  as.integer(dis),
                  beta = as.double(beta),
                                        # results -->
                  loglik = double(1), # Return value at beta
                  fail = integer(1),
                  ##DUP = TRUE,
                  PACKAGE = "eha"
                  )
        if (fit$fail) stop("Error in exp likelihood calculation")
        
        return(fit$loglik)
    }
    ncov <- NCOL(X)
    if (length(shape) == 1){
        shape <- rep(shape, ns)
    }else if (length(shape) != ns){
        stop("length(shape) must be equal to 1 or No. of strata")
    }
    
    bdim <- ns
    ncov.save <- ncov
    ncov <- 0
    beta <- numeric(bdim)
    for (i in seq_len(ns)){
        beta[i] <- log(sum(Y[, 2] - Y[, 1]) / sum(Y[, 3]))
    }

    res <- optim(beta, Fexpmin, method = "BFGS",
                 control = list(trace = as.integer(printlevel)),
                 hessian = TRUE)

    ncov <- ncov.save
    bdim <- ncov + ns
    beta <- c(rep(0, ncov), res$par)
    loglik.start <- -res$value

    res1 <- optim(beta, Fexpmin, method = "BFGS",
                 control = list(trace = as.integer(printlevel)),
                 hessian = TRUE)

    fit <- list(beta = res$par, loglik = c(loglik.start, -res$value))
    fit$fail <- res$convergence != 0
    fit$var <- try(solve(res$hessian))
    fit$shape.fixed <- TRUE
    fit$shape <- shape
    fit$shape.sd <- NULL  ## Not necessary!?!?
    ##fit$beta[bdim] <- -fit$beta[bdim] # To get "1 / lambda"! NO!!
    ##if (ncov & !center){
    
    fit$ncov <- ncov
    fit
}
