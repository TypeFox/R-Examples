##=======================================================================
## Author: Yves Deville
##
## Find ML estimate of a two-parameters 'maxlo' distribution using a
## sample 'x'. The likelihood is concentrated with respect to the shape
## parameter
##
## This distribution is a special case of a beta distribution on an
## interval (0, beta)  where 'beta' is considered as an unknown parameter.
##
## It is a reparametrisation of a GPD with location mu = 0 and negative
## shape xi < 0
##
##=======================================================================

`fmaxlo` <- function(x,
                     shapeMin = 1.25,
                     info.observed = FALSE,
                     plot = FALSE,
                     scaleData = TRUE,
                     cov = TRUE) {
    
    Cvg <- TRUE
    if (any(x <= 0)) stop("all elements in 'x' must be > 0")
    if (shapeMin <= 0) stop("'shapeMin' must be positive")
    if (shapeMin < 1) warning("'shapeMin <= 1': non-identifiable model")
    
    parnames <- c("shape", "scale")
    n <- length(x)
    M1 <- mean(x)
    
    if (scaleData) {
        ## xUnscaled <- x    ## useful ???
        x <- x / M1
        CV <- sqrt(1 - 1/n) * sd(x)
        cLogLik <- - n * log(M1)
        trans <- diag(c(1, 1 / M1))
        colnames(trans) <- rownames(trans) <- parnames
    } else {
        CV <- sqrt(1 - 1/n) * sd(x) / M1
    }
    
    if (CV > 1.00) stop("CV > 1. Estimation impossible for \"maxlo\"")
    ## if (CV > 0.999) warning("large CV value: estimation may not converge")

    ## compute an upper bound for beta. M1 <- mean(x) was done before
    M2 <- mean(x^2)
    M3 <- mean(x^3)
    
    if (scaleData) {
        betaRoots <- polyroot(c(-12*M3, 10*M3 - 6*M2, 3*(M2 - 2)))
    } else {
        betaRoots <- polyroot(c(-12*M1*M3, 10*M3 - 6*M1*M2, 3*(M2 - 2*M1^2)))
    }
    
    betaUpper <- max(Re(betaRoots))
    mind <- max(x)
    if (betaUpper < 2 * mind) betaUpper <- 2 * mind
    
    interv <- c(mind * (1 + 1e-6), betaUpper)

    ##=======================================================================
    ## o 'logLc' will be maximised. THis is the logLik concentrated with the
    ##    constraint shape >= shapeMin, sot its derivative is no longer
    ##    given by the formula for unconstrained maximisation!
    ## 
    ## o 'dlogLc' WAS used to find a zero and is used to check the
    ##    bounds in 'interval'
    ## 
    ## o other functions are needed to compute the covariance (if needed)
    ##=======================================================================
    
    ## dlogLc <- function (beta) {
    ##     xmod <- x / beta
    ##     R <- -mean(log(1.0 - xmod))
    ##     dR <- -mean( xmod / (1.0 - xmod) ) / beta
    ##     -n * ( (1.0 - R) * dR / R + 1.0 / beta ) 
    ## }
    
    logLc <- function (beta) {
        xmod <- x / beta
        R <- -mean(log(1.0 - xmod))
        rho <- 1 / R
        if (rho < shapeMin) {
            res <-  n * ( log(shapeMin) - log(beta) - (shapeMin - 1) * R) 
        } else {
            res <-  -n * ( log(R) + log(beta) - R + 1.0)
        }
        res
    }

    if (plot) cov <- TRUE
        
    if (cov) {
        
        logL <- function (parm) {
            rho <- parm[1]
            beta <- parm[2]
            xmod <- x / beta
            R <- -mean(log(1.0 - xmod))
            n*( log(rho) - log(beta) -(rho - 1.0) * R )
        }
        
        ## 'beta' is scaled, as is 'x"
        log2L <- function (alpha, beta) {
            xmod <- x / beta
            xmod1 <- 1 - xmod
            RL <- mean(-log(xmod1))
            R1 <- mean(1 / xmod1)
            R2 <- mean(1 / xmod1 / xmod1)
            
            alpha1 <- alpha - 1.0
            
            logL <- n * (log(alpha / beta)  - alpha1 * RL) 
            
            dlogL <- c(shape = n * (1 / alpha - RL),
                       scale = n * ( -1 + alpha1 * (R1 - 1) ) / beta)
            
            d2logL <- array(0, dim = c(2, 2), dimnames = list(parnames, parnames))
            d2logL["shape", "shape"] <- -n / alpha / alpha
            d2logL["shape", "scale"] <- n * (R1 - 1) / beta 
            d2logL["scale", "shape"] <- d2logL["shape", "scale"]
            d2logL["scale", "scale"] <-
                n * (alpha - alpha1 * R2 ) / beta / beta 
            
            if (scaleData) {
                logL <- logL + cLogLik
                dlogL <- trans %*% dlogL 
                d2logL <- trans %*% d2logL %*% trans
            } 
            
            list(logL = logL,
                 dlogL = dlogL,
                 d2logL = d2logL)    
        }
        
    } 
 
    ## ## check the bounds
    ## checks <- unlist(sapply(interv, dlogLc))
    ##
    ## if ( (checks[1] < 0) || (checks[2] > 0) ) {
    ##     cat("fmaxlo bounds\n"); print(checks)
    ##     warning("no interval found to maximise loglik")
    ##     Cvg <- FALSE
    ## }

    res <- optimize(f = logLc, interval = interv, maximum = TRUE)
    
    ## caution beta.hatS is for salaed data, and unscaled 'x' is lost 
    beta.hatS <- res$maximum
    alpha.hat <- -1 / mean(log(1 - x / beta.hatS))
    constr <- FALSE
    if (alpha.hat <= shapeMin) {
        constr <- TRUE
        alpha.hat <- shapeMin
    }
    
    if (scaleData) beta.hat <- M1 * beta.hatS
    else beta.hat <- beta.hatS
        
    if (!cov) {
        loglik <- res$objective 
        if (scaleData) loglik <- loglik + cLogLik
        return(list(estimate = c(shape = alpha.hat, scale = beta.hat),
                    CV = CV,
                    loglik = loglik,
                    cvg = Cvg))
    }
    
    res2 <- log2L(alpha = alpha.hat, beta = beta.hatS)
    
    if (alpha.hat > 2 && !constr) {
        
        if (info.observed) {
            info <- -res2$d2logL
            vcov <- solve(info)
            if (!inherits(cov, "try-error")) {
                vcov <- NULL
                sds <- NULL
                warning("hessian could not be inverted")
            } else {
                sds <- sqrt(diag(vcov))
            }
        } else {
            a1 <- alpha.hat - 1
            a2 <- alpha.hat - 2
            if (alpha.hat <= shapeMin + 1e-4) {
                warning("'shape' is near the bound 'shapeMin'. ",
                        "Derivatives may be missliding")
            }
            i11 <- 1 / alpha.hat / alpha.hat
            i12 <- -1 / beta.hat / a1
            i22 <- alpha.hat / a2 / beta.hat / beta.hat
            info <- matrix(n * c(i11, i12, i12, i22), nrow = 2L, ncol = 2L)
            colnames(info) <- rownames(info) <- parnames
            c11 <- alpha.hat^2 * a1^2
            c12 <- alpha.hat * a1 * a2 * beta.hat
            c22 <- a1^2 * a2 * beta.hat^2 / alpha.hat
            vcov <- matrix(c(c11, c12, c12, c22) / n, nrow = 2L, ncol = 2L)
            colnames(vcov) <- rownames(vcov) <- parnames
            sds <- sqrt(diag(vcov))
        }
        
    } else {
        ## return matrix or vector of NA
        if (constr){
            warning("'shape' is at the bound given in 'shapeMin'. ML inference results not suitable")
         } else {
             warning("'shape' is <= 2. ML inference results not suitable")
         }
        info <- matrix(NA, nrow = 2L, ncol = 2L)
        colnames(info) <- rownames(info) <- parnames
        ## print(cbind(info, -res2$d2logL)) ## for checks
        vcov <- info
        sds <- rep(NA, 2L)
        names(sds) <- parnames
    }
    
    if (plot) {
        
        dlogLc <- function (beta) {
            xmod <- x / beta
            R <- -mean(log(1.0 - xmod))
            dR <- -mean( xmod / (1.0 - xmod) ) / beta
            -n * ( (1.0 - R) * dR / R + 1.0 / beta ) 
        }
        
        if (scaleData) beta.sol <- beta.hat / M1
        else beta.sol <- beta.hat
        
        ## prepare grid and compute the limit
        lcInf <- -n * (1 + log(mean(x)))
        betas <- seq(from = interv[1], to = interv[2], length = 200)   
        fs <- sapply(betas, logLc)
        dfs <- sapply(betas, dlogLc)
        
        ind <- 1L:length(betas)
        ind <- dfs < 20
        Stext <- ifelse(scaleData, "(scaled data)", "")
        
        ## now plot logLik derivative and logLik
        opar <- par(mfrow = c(2L, 1L))
        par(mar = c(0, 5, 5, 5))
        
        ## First plot: logLik derivative. 
        plot(betas[ind], dfs[ind],
             type = "n",
             main = sprintf("'Maxlo' concentrated log-lik CV = %4.2f %s", CV, Stext),
             xlab = " ", ylab = "dlogL UNCONSTR.",
             xaxt = "n", yaxt = "n",
             xlim = interv)
        axis(side = 4)
        abline(h = 0)
        abline(v = beta.sol, col = "orangered")
        abline(v = interv, col = "darkcyan", lwd = 2)
        abline(h = 0, col = "gray")
        lines(betas[ind], dfs[ind],
              type = "l", lty = "solid", col = "red3")
    
        par(mar = c(5, 5, 0, 5))
        
        ## Second plot: logLik function 
        plot(betas[ind], fs[ind], type = "l",
             lty = "solid", col = "red3",
             xlab = "beta (scale param.)", ylab = "logL",
             xlim = interv, ylim = range(fs[ind], lcInf))
        abline(h = lcInf, col = "orchid")
        mtext(text = "lim.", side = 4, at = lcInf,
              col = "orchid")
        abline(v = interv, col = "darkcyan", lwd = 2)
        abline(v = beta.sol, col = "orangered")
        mtext(text = "betaHat", col = "orangered",
              side = 1, at = beta.sol, line = 0.5)

        par(opar)
        
    }
    
    list(estimate = c(shape = alpha.hat, scale = beta.hat),
         sd = sds,
         CV = CV,
         loglik = res2$logL,
         dloglik = res2$dlogL,
         cov = vcov,
         info = info,
         cvg = Cvg)
    
}

##=======================================================================
## Author: Yves Deville
##
## Find ML estimate of a two parameter 'maxlo' distribution WITH KNOWN
## SHAPE.
##
## NB. When the SCALE 'beta' is known, -log(1 - x / beta) is exponential
## with rate 'alpha'.
##
## TODO: add the 'scaleData' argument, find bounds for the optim...
##
##=======================================================================

`fmaxlo1` <- function(x,
                      shape = 1.5,
                      plot = FALSE) {
  
  if (any(x <= 0)) stop("all elements in 'x' must be >0")  
  if (shape < 1) stop("'shape <= 1': non-identifiable model")
  
  n <- length(x)
  parnames <- c("shape", "scale")

  CV <- sqrt(1 - 1/n) * sd(x) / mean(x) 
  if (CV > 0.99) warning("large CV value: estimation may not converge")

  alpha <- shape
  alpha1 <- alpha - 1
  
  ## could be used to find a zero
  dlogL1 <- function (beta) {
    xmod <- x / beta
    R <- -mean(log(1.0 - xmod))
    S1 <- mean( xmod / (1.0 - xmod) )
    n * ( -1  + alpha1 * S1 ) /beta 
  }
  
  logL1 <- function (beta) {
    xmod <- x/beta
    R <- -mean(log(1.0 - xmod))
    n * (log(alpha / beta) - alpha1 * R)
  }
  
  mind <- max(x)
  interv <- c(mind + 1e-6, 10 * mean(x))
  checks <- unlist(sapply(interv, dlogL1))
  if ( (checks[1] < 0) ||  (checks[2] > 0) ) 
    warning("no interval found to maximise loglik")
  res <- optimize(f = logL1, interval = interv, maximum = TRUE)
  beta.hat <- res$maximum
  loglik <- res$objective
  dloglik <- dlogL1(beta.hat)
  
  if (alpha >= 2) {
    info <-  n * alpha / (alpha-2) / beta.hat / beta.hat
    cov <- solve(info)
    sds <- sqrt(diag(cov))
  } else {
    warning("'shape' is < 2 ML inference results not suitable")
    info <- NULL
    cov <- NULL
    sds <- NULL
  }
  
  if (plot) {
    betas <- seq(from = mind + 1e-5, to = 10*mean(x), length = 200)
    fs <- sapply(betas, logL1)
    dfs <- c(NA, diff(fs) / diff(betas))
    ind <- 1L:length(betas)
    ind <- dfs < 20
    plot(betas[ind], dfs[ind], type = "l", lty = "dotted")
    lines(betas[ind], sapply(betas[ind], dlogL1), col = "red")
    abline(v = mind, col = "purple")
    abline(h = 0)
    abline(v = beta.hat)
  }
  
  list(estimate = c(scale = beta.hat),
       sd = sds,
       CV = CV,
       loglik = loglik,
       dloglik = dloglik,
       cov = cov,
       info = info)
  
}

