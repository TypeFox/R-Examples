glmmbootFit <- function (X, Y, weights = rep(1, NROW(Y)), 
                         start.coef = NULL,
                         cluster = rep(1, length(Y)),                        
                         offset = rep(0, length(Y)),
                         family = binomial(),
                         control = list(epsilon = 1.e-8, maxit = 200,
                           trace = FALSE), 
                         boot = 0){

    if (is.list(control)) {
        if (is.null(control$epsilon))
          control$epsilon <- 1e-08
        if (is.null(control$maxit))
          control$maxit <- 200
        if (is.null(control$trace))
          control$trace <- FALSE
    }
    else {
        stop("control must be a list")
    }

    X <- as.matrix(X)

    nobs <- NROW(Y)
    if (family$family == "binomial"){ # This will be better later!
        ## From 'binomial':
        if (NCOL(Y) == 1) {
            if (is.factor(Y)) Y <- Y != levels(Y)[1]
            n <- rep.int(1, nobs)
            if (any(Y < 0 | Y > 1)) stop("Y values must be 0 <= Y <= 1")
            ##mustart <- (weights * Y + 0.5)/(weights + 1)
            m <- weights * Y
            if (any(abs(m - round(m)) > 0.001))
              warning("non-integer #successes in a binomial glm!")
        } else if (NCOL(Y) == 2) {
            if (any(abs(Y - round(Y)) > 0.001))
              warning("non-integer counts in a binomial glm!")
            n <- Y[, 1] + Y[, 2]
            Y <- ifelse(n == 0, 0, Y[, 1]/n)
            weights <- weights * n
            ##mustart <- (n * Y + 0.5)/(n + 1)
        } else
        stop("for the binomial family, Y must be a vector of 0 and 1's\n",
             "or a 2 column matrix where col 1 is no. successes and col 2 is no. failures")
        ## End of 'from binomial'.
    }else{
        if (NCOL(Y) != 1)
          stop("It is only the binomial family that allows two columns for Y")
    }
    coli <- match("(Intercept)", colnames(X))
    with.intercept <- !is.na(coli)

    if (is.null(offset)) offset <- rep(0, length(Y))
    glmFit <- glm.fit(X, Y, weights = weights,
                      ##start = c(0, start.coef),
                      offset = offset,
                      family = family,
                      control = control,
                      intercept = with.intercept,
                      )
    predicted <- glmFit$fitted.values
    cluster.null.deviance <- glmFit$deviance
    if (with.intercept)
      X <- X[, -coli, drop = FALSE]
    p <- ncol(X)
    if (is.null(start.coef)){
        start.coef <- numeric(p) # Start values equal to zero,
        if (FALSE){## Not sensible below:
            if (family$family == "binomial"){
                start.coef[1] <- log(mean(Y) / (1 - mean(Y)))
            }else if (family$family == "poisson"){
                start.coef[1] <- log(mean(Y))
            }else{ ## this is a proviso!!
                start.coef[1] <- mean(Y)
            }
        } ## End 'if (FALSE)'
    }else{                   
        if (length(start.coef) != p) stop("beta.start has wrong length")
    }

    ord <- order(cluster)

    Y <- Y[ord]
    X <- X[ord, ,drop = FALSE]
    cluster <- cluster[ord]

    if (family$family == "binomial"){
        if (family$link == "logit"){
            fam <- 0
        }else if (family$link == "cloglog"){
            fam <- 1
        }else{
            stop("Unknown link function; only 'logit' and 'cloglog' implemented")
        }
    }else if (family$family == "poisson"){
        fam <- 2
    }else{
        stop("Unknown family; only 'binomial' and 'poisson' implemented")
    }
    
    famSize <- as.vector(table(cluster))
    nFam <- length(famSize)
  
    ## cat("nFam = ", nFam, "\n")

    if (p >= 1){
        means <- colMeans(X)
        X <- scale(X, center = TRUE, scale = FALSE)
        ## cat("means = ", means, "\n")
        fit <- .C("glmm_boot",
                  as.integer(fam),
                  as.integer(p),
                  as.double(start.coef),
                  as.integer(cluster),
                  as.double(weights),
                  as.double(t(X)),       # Note! #
                  as.double(Y),
                  as.double(offset),
                  as.integer(famSize),
                  as.integer(nFam),
                  as.double(control$epsilon),
                  as.integer(control$maxit),
                  as.integer(control$trace),
                  as.integer(boot),
                  beta = double(p),
                  predicted = as.double(predicted), # Watch up! #
                  fitted = double(length(Y)),
                  loglik = double(1),
                  variance = double(p * p),
                  info = integer(1),
                  frail = double(nFam),
                  bootP = double(1),
                  bootLog = double(boot),
                  convergence = integer(1),
                  ## DUP = FALSE,
                  PACKAGE = "eha"
                  )
        fit$frail[fit$frail < -999] <- -Inf
        fit$frail[fit$frail > 999] <- Inf
        res <- list(coefficients = fit$beta,
                    predicted = fit$predicted,
                    fitted = fit$fitted,
                    logLik = fit$loglik,
                    cluster.null.deviance = cluster.null.deviance,
                    ## Corrects for "centering" above (Added 0.72):
                    ## frail = fit$frail - means * fit$beta,
                    ## Should be (The above is obviously wrong):
                    frail = fit$frail - sum(means * fit$beta),
                    ## ?? Changed 2011-08-08. 
                    bootLog = fit$bootLog,
                    bootP = fit$bootP,
                    info = fit$info)
        if (!fit$info){
            res$variance <- matrix(fit$variance, nrow = p, ncol = p)
            ## Removed minus sign above, 0.65-4
            res$sd <- sqrt(diag(res$variance))
        }else{
            res$variance <- NULL
            res$sd <- rep(NA, p)
        }
        res$boot_rep <- boot
        
        return(res)
    }else{ # A null model:
        
        fit <- .C("glmm_boot0",
                  as.integer(fam),
                  ##as.integer(p),
                  ##as.double(start.coef),
                  as.integer(cluster),
                  as.double(weights),       ## Note! ##
                  as.double(Y),
                  as.double(offset),
                  as.integer(famSize),
                  as.integer(nFam),
                  ##as.double(control$epsilon),
                  ##as.integer(control$maxit),
                  as.integer(control$trace),
                  as.integer(boot),
                  predicted = as.double(predicted),
                  ##beta = double(p),
                  fitted = double(length(Y)),
                  loglik = double(1),
                  ##hessian = double(p * p),
                  frail = double(nFam),
                  bootP = double(1),
                  bootLog = double(boot),
                  convergence = integer(1),
                  ## DUP = FALSE,
                  PACKAGE = "eha"
                  )
        res <- list(coefficients = NULL,
                    predicted = fit$predicted,
                    fitted = fit$fitted,
                    logLik = fit$loglik,
                    cluster.null.deviance = cluster.null.deviance,
                    frail = fit$frail,
                    bootLog = fit$bootLog,
                    bootP = fit$bootP)
        res$variance <- NULL
        res$sd <- NULL
        res$boot_rep <- boot
        
        return(res)
    }
    
}
