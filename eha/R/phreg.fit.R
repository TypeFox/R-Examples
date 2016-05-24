phreg.fit <- function(X,
                      Y,
                      dist,
                      strata,
                      offset,
                      init,
                      shape,
                      control,
                      center = NULL){# This 'center' not used.

    if (dist == "weibull"){
        dis <- 0
        center <- TRUE
        intercept <- FALSE
    }else if(dist == "loglogistic"){
        dis <- 1
        center <- FALSE
        intercept <- TRUE
        ##namn <- colnames(X)
        ##X <- cbind(rep(1, NROW(X)), X)
        ##colnames(X) <- c("(Intercept)", namn)
        if (!missing(init)) init <- c(0, init)
    }else if (dist == "lognormal"){
        dis <- 2
        center <- FALSE
        intercept <- TRUE
        ##namn <- colnames(X)
        ##X <- cbind(rep(1, NROW(X)), X)
        ##colnames(X) <- c("(Intercept)", namn) 
        if (!missing(init)) init <- c(0, init)
    }else if (dist == "ev"){
        dis <- 3
        center <- TRUE
        intercept <- FALSE
    }else if (dist == "gompertz"){
        stop("phreg.fit cannot be used with 'gompertz', try 'gompreg'")
    }else if (dist == "pch"){
        stop("phreg.fit cannot be used with 'pch', try 'pchreg'")
    }else{
        stop(paste(dist, "is not an implemented distribution"))
    }

    nn <- NROW(X)
    ncov <- NCOL(X)

    
    ## Should we really have _weighted_ means? Yes! (Cf. coxreg.fit)
    if (ncov && center){
        wts <- Y[, 2] - Y[, 1]
        w.means <- apply(X, 2, weighted.mean, w = wts)
        means <- colMeans(X)
        ##means <- apply(X, 2, mean)
        for (i in 1:ncov){
            X[, i] <- X[, i] - means[i]
        }
    }
    if (missing(strata) || is.null(strata)){
        strata <- rep(1, nn)
        ns <- 1
    }else{
        strata <- as.integer(factor(strata))
        ns <- max(strata)
    }

    if (length(strata) != nn) stop("Error in stratum variable")
    if (missing(offset) || is.null(offset))
        offset <- rep(0, nn)

    if (missing(init) || is.null(init))
        init <- rep(0, ncov)
    if (length(init) != ncov) stop("Error in init")

    printlevel <- control$trace
    iter <- control$maxiter


    nstra <- c(0, cumsum(table(strata)))
    if (all(shape <= 0)){ ## Then shape is estimated in all strata

        bdim <- ncov + 2 * ns
        if (ns > 0){
            ord <- order(strata)
            X <- X[ord, , drop = FALSE]
            Y <- Y[ord, , drop = FALSE]
            offset <- offset[ord]
        }

            fit <- .C("phsup",
                  iter = as.integer(iter), #maxit on ip, actual on op.
                  as.double(control$eps),
                  as.integer(printlevel),
                                        #
                  as.integer(ns), # No. of strata
                  as.integer(nstra),
                  as.integer(nn),
                  as.integer(ncov),
                  as.integer(bdim),
                                        #
                  as.double(Y[, 1]),  ## 'enter'
                  as.double(Y[, 2]),  ## 'exit'
                  as.integer(Y[, 3]), ## 'event'
                                        #
                  ##as.double(t(scale(X, center = TRUE, scale = FALSE))),
                  as.double(t(X)), # NOTE; scaling already done!
                  ## NOTE transpose!
                  as.double(offset),
                                        #
                  as.integer(dis),     # baseline distribution
                  as.double(init),     # 'start.beta'
                  beta = double(bdim), # results -->
                  lambda = double(ns),
                  lambda.sd = double(ns),
                  shape = double(ns),  ## "p"
                  shape.sd = double(ns),
                                        #
                  loglik = double(2), # [1] == start, [2] == maximized
                  dloglik = double(bdim),
                  variance = double(bdim * bdim),
                  sctest = double(1),
                                        #
                  conver = integer(1),
                  fail = integer(1),
                  ## DUP = FALSE,
                  PACKAGE = "eha"
                  )

        if (fit$fail) return(list(fail = fit$fail,
                                  n.strata = ns,
                                  value = fit$beta[fit$fail])
                             )
        ##if (!center){ # Transform to original values
        ##if (!center){ # NOPE!!! Transform to original values
        if (ncov && center){ # 'center' is deprecated' ... (not reported)
            dxy <- diag(2 * ns + ncov)
            for (i in 1:ns){ ## Really a HACK ??!!!!!!!!!!!!!!!
                row <- ncov + 2 * i - 1
                col <- row + 1
                ## fit$beta[row] <- -fit$beta[row] NOT ANY MORE!
                if (ncov){
                    pi.hat <- exp(fit$beta[col])
                    scale.corr <- sum(means * fit$beta[1:ncov]) /
                        pi.hat
                    fit$beta[row] <- fit$beta[row] + scale.corr
                    dxy[row, 1:ncov] <- means / pi.hat
                    dxy[row, col] <- -scale.corr
                }
                ##dxy[row, row] <- -1
            }
        }
        coef.names <- colnames(X)
        if (ns > 1){
            for (i in 1:ns){
                coef.names <- c(coef.names,
                                paste("log(scale)", as.character(i), sep =":"),
                                paste("log(shape)", as.character(i), sep =":"))
            }

        }else{
            coef.names <- c(coef.names,
                            "log(scale)", "log(shape)")
        }


        fit$shape.fixed <- FALSE

    }else{  ## Then shape is fixed in all strata:

        if (ns >= 2) warning("'strata' is maybe useful for 'fixed shape' regression.\n Maybe include stratum variable as a factor in the model instead.")
        bdim <- ncov + ns

        if (length(shape) == 1){
            shape <- rep(shape, ns)
        }else if (length(shape) != ns){
          stop("length(shape) must be equal to 1 or No. of strata")
        }

        fit <- .C("phexpsup",
                  iter = as.integer(iter), #maxit on ip, actual on op.
                  as.double(control$eps),
                  as.integer(printlevel),
                                        #
                  as.integer(ns), # No. of strata
                  as.integer(nstra),
                  as.integer(nn),
                  as.integer(ncov),
                  as.integer(bdim),
                                        #
                  as.double(Y[, 1]),  ## 'enter'
                  as.double(Y[, 2]),  ## 'exit'
                  as.integer(Y[, 3]), ## 'event'
                                        #
                  ##as.double(t(scale(X, center = TRUE, scale = FALSE))),
                  as.double(t((X))), #NOTE: scaling already done!
                  as.double(offset),
                  as.double(shape), ## "p" (fixed)
                  as.integer(dis),

                  as.double(init),     # 'start.beta'
                  beta = double(bdim), # results -->
                  lambda = double(ns),
                  lambda.sd = double(ns),
                                        #
                  loglik = double(2), # [1] == start, [2] == maximized
                  dloglik = double(bdim),
                  variance = double(bdim * bdim),
                  sctest = double(1),
                                        #
                  conver = integer(1),
                  fail = integer(1),
                  ## DUP = FALSE,
                  PACKAGE = "eha"
                  )
        if (fit$fail) return(list(fail = fit$fail,
                                  n.strata = ns,
                                  value = fit$beta[fit$fail])
                             )
        fit$shape.fixed <- TRUE
        fit$shape <- shape
        fit$shape.sd <- NULL  ## Not necessary!?!?
        ##fit$beta[bdim] <- -fit$beta[bdim] # To get "1 / lambda"! NO!!
        ##if (!center){# i.e., never...
        if (center && ncov){# always transform back ....
            dxy <- diag(bdim)
            if (ncov){
                dxy[bdim, 1:ncov] <- means / shape
                scale.corr <- sum(means * fit$beta[1:ncov]) / shape
                fit$beta[bdim] <-
                    fit$beta[bdim] + scale.corr
            }
            ##dxy[bdim, bdim] <- -1
        }
        coef.names <- c(colnames(X), "log(scale)")

        ## Note; this is really a "hack"!!!!!!!!!!!!!!!
    }

    ##cat("done!\n")

    if (!fit$fail){
        if (center && ncov){ # Transform back...
            var <- dxy %*% matrix(fit$variance, bdim, bdim) %*% t(dxy)
        }else{
            var <- matrix(fit$variance, bdim, bdim)
       }
        colnames(var) <- rownames(var) <- coef.names
    }
    else
        var <- NULL

    coefficients <- fit$beta
    names(coefficients) <- coef.names
    if (intercept) df <- ncov - 1 else df <- ncov
    list(coefficients = coefficients,
         df = df,
         var = var,
         loglik = fit$loglik,
         score = fit$sctest,
         conver = fit$conver,
         fail = fit$fail,
         iter = fit$iter,
         n.strata = ns,
         shape = fit$shape
         )

}
