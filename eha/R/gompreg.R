gompreg <- function(X, Y, strata, offset, init, control, center){
    ## Gompertz proportional hazards:
    ## Stratum i: h_i(t; beta) = a_i * exp(t / b_i) * exp(x*beta),
    ## i = 1, ..., ns
    ## NOTE: Y is nn x 3!

    ## This is the canonical version.
    
    if (!is.matrix(X)) X <- matrix(X, ncol = 1)
    nn <- NROW(X)
    ncov <- NCOL(X)
    if (NROW(Y) != nn) stop("Y NROW error")
    if (NCOL(Y) != 3) stop("Y NCOL error", )
    if (ncov){
        wts <- Y[, 2] - Y[, 1]
        means <- apply(X, 2, weighted.mean, w = wts)
        means <- apply(X, 2, mean) ## ??

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

    bdim <- ncov + 2 * ns

    dGomp <- function(beta){
        ## Calculates the first derivatives of a Gompertz regression (stratified)
        ## beta[1:ncov] = the beta coefficients.
        ## beta[ncov + 1], beta[ncov + 3], ... = scale[1, 2, ...] ("gamma")
        ## beta[ncov + 2], beta[ncov + 4], ... = shape[1, 2, ...] ("alpha")
        
        enter <- Y[, 1]
        exit <- Y[, 2]
        event <- Y[, 3]
        if (ncov){
            bz <- offset + X %*% beta[1:ncov]
        }else{
            bz <- offset
        }
        
        grad <- numeric(ncov + 2 * ns)
        
        for (i in 1:ns){
            T0 <- enter[strata == i]
            T <- exit[strata == i]
            D <- event[strata == i]
            z <- X[strata == i, ,drop = FALSE]
            
            ezb <- exp(bz[strata == i])
            alpha <- beta[ncov + 2 * i] # shape
            gamma <- beta[ncov + 2 * i - 1] # scale
            scale <- exp(gamma)
            shape <- exp(alpha)
            emgamma <- exp(-gamma)
            
            ezbH <- ezb * 
                (Hgompertz(T, scale = scale, shape = shape,
                           param = "canonical") -
                 Hgompertz(T0, scale = scale, shape = shape,
                           param = "canonical"))
            
            grad[ncov + 2 * i] <- sum(D) -  sum(ezbH)
            
            grad[ncov + 2 * i - 1] <- -sum(D * T) * emgamma - sum(D) +
                sum(ezb *
                    (T * hgompertz(T, scale = scale, shape = shape,
                                   param = "canonical") -
                     T0 * hgompertz(T0, scale = scale, shape = shape,
                                    param = "canonical")))
            
            if (ncov){
                for (j in 1:ncov){
                    grad[j] <- grad[j] + sum(D * z[, j]) - sum(ezbH * z[, j])
                }
            }
        }
        grad
    }    
    
    Fmin <- function(beta){
        
        total <- 0
        
        for (i in 1:ns){
            scale <- exp(beta[ncov + 2 * i - 1])
            shape <- exp(beta[ncov + 2 * i]) ## Added!
            if (ncov){
                bz <- offset[strata == i] +
                    X[strata == i, , drop = FALSE] %*% beta[1:ncov]
            }else{
                bz <- offset[strata == i]
            }
            ebz <- exp(bz)
            H1 <- Hgompertz(Y[strata == i, 2], shape = shape, scale = scale,
                            param = "canonical")
            H0 <- Hgompertz(Y[strata == i, 1], shape = shape, scale = scale,
                            param = "canonical")
            h <- hgompertz(Y[strata == i, 2], shape = shape, scale = scale,
                           log = TRUE, param = "canonical")
            ret1 <- sum(Y[strata == i, 3] * (h + bz))
            ret2 <- sum(ebz * (H1 - H0))
            total <- total + ret1 - ret2
        }

        return(total)
    }

    
    ## First, fit the 'null' model:
    beta0 <- numeric(2 * ns)
    ncov.save <- ncov
    ncov <- 0
    ## Start values for null model:
    for (i in 1:ns){
        enter <- Y[strata == i, 1]
        exit <- Y[strata == i, 2]
        event <- Y[strata == i, 3]
        score <- offset[strata == i]
        beta0[(2 * i - 1):(2 * i)] <-
            gompstartCanonical(enter, exit, event, score)
        ##cat("\nSTART values (scale, shape):",
          ##  beta[(ncov + 2 * i - 1):(ncov + 2 * i)], "\n\n")
    }
    res0 <- optim(beta0, Fmin, gr = dGomp,
                 method = "BFGS",
                 control = list(fnscale = -1, reltol = 1e-10),
                 hessian = FALSE)
    ## Done; now the real thing:
    ncov <- ncov.save
    
    beta <- numeric(bdim)
    if (ncov)
        beta[1:ncov] <- init  # Start values
    ## Start values
    for (i in 1:ns){
        enter <- Y[strata == i, 1]
        exit <- Y[strata == i, 2]
        event <- Y[strata == i, 3]
        score <- offset[strata == i] + X[strata == i, ,drop = FALSE] %*% init
        beta[(ncov + 2 * i - 1):(ncov + 2 * i)] <-
            gompstartCanonical(enter, exit, event, score)
        ##cat("\nSTART values (scale, shape):",
          ##  beta[(ncov + 2 * i - 1):(ncov + 2 * i)], "\n\n")
    }
    

    res <- optim(beta, Fmin, gr = dGomp,
                 method = "BFGS",
                 control = list(fnscale = -1, reltol = 1e-10),
                 hessian = TRUE)
    if (res$convergence != 0) stop("[gompreg]: No convergence")
    coefficients <- res$par
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

    names(coefficients) <- coef.names
    fit <- list(coefficients = coefficients,
                loglik = c(res0$value, res$value) # c(res0$value, res$value)
                )
    fit$gradient <- dGomp(coefficients)
    fit$pfixed <- FALSE
    fit$var <- tryCatch(solve(-res$hessian), error = function(e) e)

    fit$hessian <- res$hessian
    if (is.matrix(fit$hessian)){
        colnames(fit$hessian) <- coef.names
        rownames(fit$hessian) <- coef.names
    }

    ## Fixing the shape parameter(s) due to centering:
    ##cat(" Before: fit$coefficients = ", fit$coefficients, "\n")
    ##cat("means = ", means, "\n")
    ##cat("coefficients = ", fit$coefficients[1], "\n")
    if (ncov){
        dxy <- diag(2 * ns + ncov)
        for (i in 1:ns){
            row <- ncov + 2 * i
            shape.corr <- sum(means * fit$coefficients[1:ncov]) #/
            fit$coefficients[row] <- fit$coefficients[row] - shape.corr
            dxy[row, 1:ncov] <- means
        }
        if (is.numeric(fit$var)) fit$var <- dxy %*% fit$var %*% t(dxy)
        fit$dxy <- dxy
    }
    ##cat(" After: fit$coefficients = ", fit$coefficients, "\n")
    if (is.matrix(fit$var)){
        colnames(fit$var) <- coef.names
        rownames(fit$var) <- coef.names
    }
    
    fit$n.strata <- ns
    fit$df <- ncov
    fit$fail <- FALSE # Optimist!
    fit$param <- "canonical"
    fit
}
                      
