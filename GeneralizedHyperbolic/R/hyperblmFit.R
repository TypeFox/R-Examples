hyperblmFit <- function(x, y, paramStart = NULL, offset = NULL,
                        method = c("Nelder-Mead","BFGS","nlm"),
                        startMethod = c("Nelder-Mead","BFGS"),
                        startStarts = c("BN","US","FN","SL","MoM"),
                        maxiter = 100, tolerance = 0.001,
                        controlBFGS = list(maxit = 1000),
                        controlNM = list(maxit = 1000),
                        maxitNLM = 10000,
                        controlCO = list(), silent = TRUE,
                        breaks = NULL, ...)
{
    if(is.null(n <- nrow(x)))
        stop("'x' must be a matrix")
    if(n == 0)
        stop("0 (non-NA) cases")
    p <- ncol(x)
    if(p == 0){
        return(list(coefficients = numeric(0), residuals = y,
                    fitted.values = 0*y, rank = 0, df.residual = length(y)))
    }else{
        xNames <- colnames(x)
    }
    ny <- NCOL(y)
    if(is.matrix(y) && ny == 1)
        y <- drop(y)
    if(!is.null(offset))
        y <- y - offset
    if(NROW(y) != n)
        stop("incompatible dimensions")
    storage.mode(x) <- "double"
    storage.mode(y) <- "double"
    ## set default error message
    errMessage <- ""


    if(is.null(paramStart)){
        qrx <- qr(x)
        resids <- qr.resid(qrx, y)
        Beta <- as.numeric(qr.coef(qrx, y))
        startInfo <- hyperbFitStand(resids, startMethod = startMethod,
                                    method = "constrOptim",
                                    startValues = startStarts,
                                    silent = TRUE, ...)
        residsParamStart <- as.numeric(startInfo$param)

        ## change residsParamStart to param set number 1
        ## (mu,delta,pi,zeta)
        distparam <-
            as.numeric(hyperbChangePars(2, 1, param = residsParamStart))[-1]
        coef <- c(residsParamStart[1] + Beta[1], Beta[-1])
        breaks <- startInfo$breaks
    } else {
        if(length(paramStart)!= (3 + p))
            stop(paste("Parameters start value should be of dimension",
                       3 + p, sep = " "))
        if(paramStart[2] <= 0)
            stop("zeta in paramStart must be greater than zero")
        if(paramStart[1] <= 0)
            stop("delta in paramStart must be greater than zero")
        distparam <- paramStart[1:3]
        coef <- paramStart[-(1:3)]
    }

    ## Set some parameters to help with optimization
    eps <- 1e-16
    ratioCoef <- 1
    rationParam <- 1
    iter <- 0

    ## A artificial iterations counter that equals to the number of iterations
    ## when the optimization does not converge and is greater than maximum
    ## iterations number when the optimization converges.
    ## This helps to break the loop.
    hiter <- 0

    sOnellfunc <- function(coef){
        KNu <- besselK(distparam[3], nu = 1)
        resids <- y - as.vector(x %*% as.matrix(coef))
        hyperbDens <- (2*distparam[1]*sqrt(1 + distparam[2]^2)*KNu)^(-1)*
                          exp(-distparam[3]*(sqrt(1 + distparam[2]^2)*
                              sqrt(1 + ((resids)/distparam[1])^2) -
                              distparam[2] * ((resids)/distparam[1])))
        ##cat("log-likelihood is", -sum(log(hyperbDens)), "\n")
        return(-sum(log(hyperbDens)))
    }


    sTwollfunc <- function(distparam) {
        ## Protect against attempts to make parameters < 0
        if (distparam[3] <= eps) return(1e99)
        if (distparam[1] <= eps) return(1e99)
        KNu <- besselK(distparam[3], nu = 1)
        resids <- y - as.vector(x %*% as.matrix(coef))
        hyperbDens <- (2*distparam[1]*sqrt(1 + distparam[2]^2)*KNu)^(-1)*
                          exp(-distparam[3]*(sqrt(1 + distparam[2]^2)*
                              sqrt(1 + ((resids)/distparam[1])^2) -
                              distparam[2]*((resids)/distparam[1])))
        ##cat("log-likelihood is", -sum(log(hyperbDens)), "\n")
        return(-sum(log(hyperbDens)))
    }


    while ( !is.null(coef) && !is.null(distparam) && hiter < maxiter)
    {
        output <- numeric(7)
        ind <- 1:6

        ##Stage one find the coefficient
        ##optimize method
        if(!is.null(distparam)){
            coefOld <- coef
            if(method == "BFGS"){
                tryOpt <- try(optim(coef, sOnellfunc, NULL, method = "BFGS",
                                    control = controlBFGS, ...),
                              silent = silent)
                if (class(tryOpt) == "try-error"){
                    errMessage <- unclass(tryOpt)
                } else {
                    optOutCoef <- tryOpt
                }
            }
            if(method == "Nelder-Mead"){
                tryOpt <- try(optim(coef, sOnellfunc, NULL,
                                    method = "Nelder-Mead",
                                    control = controlNM, ...),
                              silent = silent)
                if (class(tryOpt) == "try-error"){
                    errMessage <- unclass(tryOpt)
                } else {
                    optOutCoef <- tryOpt
                }
            }

            if(method == "nlm"){
                ind <- c(2, 1, 5, 4)
                tryOpt <- try(nlm(sOnellfunc, coef,
                                  iterlim = maxitNLM, ...),
                              silent = silent)
                if (class(tryOpt) == "try-error"){
                    errMessage <- unclass(tryOpt)
                } else {
                    optOutCoef <- tryOpt
                }
            }
            if (errMessage == ""){
                coef <- as.numeric(optOutCoef[[ind[1]]])
                ratioCoef <- max(abs((coefOld - coef)/coefOld))
            }
            else coef <- NULL
        }

        ## Stage Two Optmization
        if (!is.null(coef)){
            ind <- 1:6
            distparamOld <- distparam
            tryOpt <- try(optOut<- constrOptim(theta = distparam,
                                               sTwollfunc, NULL,
                                               ui = diag(c(1, 0, 1)),
                                               ci = c(0, -1e+99, 0),
                                               control = controlCO, ...),
                          silent = silent)

            if (class(tryOpt) == "try-error"){
                errMessage <- unclass(tryOpt)
            } else {
                optOut <- tryOpt
            }
            if (errMessage == ""){
                distparam <- as.numeric(optOut[[ind[1]]])
                ratioParam <- max(abs((distparamOld - distparam)/distparamOld))
            } else {
                distparam <- NULL
            }
        }
        iter <- iter + 1
        if(ratioCoef <= tolerance && ratioParam <= tolerance )
            hiter <- maxiter + 1
        else hiter <- iter
    }


    if (!is.null(distparam)){
        alpha <- distparam[3] * sqrt(1 + distparam[2]^2) / distparam[1]
        beta <- distparam[3] * distparam[2] / distparam[1]
        ordistparam <- c(distparam[1], alpha, beta)

        KNu <- besselK(distparam[3], nu = 1)
        resids <- y - as.vector(x %*% as.matrix(coef))
        hyperbDens <- (2*distparam[1]*sqrt(1 + distparam[2]^2)*KNu)^(-1)*
                          exp(-distparam[3]*(sqrt(1 + distparam[2]^2)*
                              sqrt(1 + (resids/distparam[1])^2) -
                              distparam[2] * (resids/distparam[1])))

        ## Calculate the MLE
        maxLik <- sum(log(hyperbDens))
        if(ratioCoef <= tolerance && ratioParam <= tolerance)
            conv <- 0
        else if(ratioCoef > tolerance && ratioParam <= tolerance)
            conv <- 1
        else if(ratioCoef <= tolerance && ratioParam > tolerance)
            conv <- 2
        else conv <- 3
        iter <- iter
    } else {
        ordistparam <- NULL
        maxLik <- NULL
        conv <- 3
        iter <- NULL
    }


    ## Fitted value
    fits <- x %*% as.matrix(coef)

    if(!is.null(offset)){
        fits <- fits + offset
    }

    resids <- y - fits
    m.r <- mean(resids)

    fits <- fits + m.r
    resids <- resids - m.r

    distributionParams <- c(-m.r, ordistparam)
    coef[1] <- coef[1] + m.r

    regressionResult <- list(coefficients = coef,
                             distributionParams = distributionParams,
                             MLE = maxLik, method = method,
                             convergence = conv, iterations = iter,
                             fitted.values=fits,
                             paramStart = paramStart, breaks = breaks,
                             residsParamStart = residsParamStart,
                             xNames = xNames,
                             residuals = resids,
                             xMatrix = x, yVec = y)

    regressionResult
}
