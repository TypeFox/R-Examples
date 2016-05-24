# some simple R wrapper functions

crossprodF <- function(X, Y, ...)
{
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  out <- .Fortran("crossprod",
                  X = X,
                  Y = Y,
                  n = as.integer(nrow(X)),
                  p = as.integer(ncol(X)),
                  q = as.integer(ncol(Y)),
                  XTY = matrix(0, ncol(X), ncol(Y)) )
  return(out$XTY)
}


cov.wtF <- function(X, Z, ...)
{
    X <- as.matrix(X)
    Z <- as.matrix(Z)
    n <- nrow(X)
    p <- ncol(X)
    nz <- nrow(Z)
    G <- ncol(Z)
    if ( n != nz ) stop("X and Z must have same number of rows")

    tmp <- .Fortran("covwt",
                    x = as.double(X),
                    z = as.double(Z),
                    n = as.integer(n),
                    p = as.integer(p),
                    G = as.integer(G),
                    mu = double(p*G),
                    W = double(p*p*G) )
    out <- list( mu = matrix(tmp$mu, p,G), W = array(tmp$W, c(p,p,G)) )
    return(out)
}


##############################################################################
###                               EVV model                               ####
##############################################################################
emEVV <- function(data, parameters, prior = NULL, control = emControl(),
                  warn = NULL, ...)
{
    z <- estepEVV(data, parameters = parameters, warn = warn)$z
    meEVV(data, z = z, prior = prior, control = control,
          Vinv = parameters$Vinv, warn = warn)
}

####
meEVV <- function(data, z, prior = NULL, control = emControl(),
                  Vinv = NULL, warn = NULL, ...)
{
    if(is.null(warn)) warn <- mclust.options("warn")
    dimdat <- dim(data)
    oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
    if(oneD || length(dimdat) != 2)
        stop("data should in the form of a matrix")
    data <- as.matrix(data)
    n <- nrow(data)
    p <- ncol(data)
    z <- as.matrix(z)
    dimz <- dim(z)
    if(dimz[1] != n)
        stop("data and z should have the same row dimension")
    K <- dimz[2]
    if (!is.null(Vinv)) {
        G <- K - 1
        if(Vinv <= 0) Vinv <- hypvol(data, reciprocal = TRUE)
    } else G <- K
    #
    if(all(is.na(z))) {
        WARNING <- "z is missing"
        if(warn) warning(WARNING)
        variance <- list(modelName = "EVV", d = p, G = G,
                         scale = NA, shape = rep(NA,p), orientation = array(NA,c(p,p,G)))
        parameters <- list(pro=rep(NA,G), mean=matrix(NA,p,G), 
                           variance=variance, Vinv=Vinv)
        return(structure(list(modelName="EVV", prior=prior, n=n, d=p,
                              G=G, z=z, parameters=parameters,
                              control=control, loglik=NA),
                         WARNING = WARNING, returnCode = 9))
        
    }
    if(any(is.na(z)) || any(z < 0) || any(z > 1))
        stop("improper specification of z")
    lwork <- max(3 * min(n, p) + max(n, p), 5 * min(n, p))
    storage.mode(z) <- "double"
    #
    #
    # MICHAEL from here------------------------------------------
    #
    # without prior specification
    if(is.null(prior)) {
        temp <- .Fortran( "meevv",
                          x = as.double(data),
                          z = as.double(z),
                          n = as.integer(n),
                          p = as.integer(p),
                          G = as.integer(G),
                          Gnoise = as.integer(K),
                          mu = double(p*G),
                          O = double(p*p*G),
                          shape.o = double(p*p*G),
                          scale = double(G),
                          shape = double(p*G),
                          pro = double(K),
                          Vinv = as.double( if (is.null(Vinv)) -1 else Vinv ), 
                          loglik = double(1),
                          eqpro = as.logical(control$equalPro),
                          itmax = as.integer(control$itmax[1]),
                          tol = as.double(control$tol[1]),
                          eps = as.double(control$eps),
                          niterout = integer(1),
                          errout = double(1),
                          lwork = as.integer(lwork),
                          info = as.integer(0),
                          package = "mclust")
    } else {
        # with prior
        priorParams <- do.call(prior$functionName,
                               c(list(data = data, G = G, modelName = "EVV"),
                                 prior[names(prior) != "functionName"]))
        # temp <- .Fortran("meevvp", ...)
        temp <- list(x = data, z = z, n = n, p = p, G = G,
                     mu = double(p*G), O = double(p*p*G), shape.o = double(p*p*G),
                     scale = double(G), shape = double(p*G), pro = double(K),
                     loglik = NA, eqpro = as.logical(control$equalPro),
                     itmax = as.integer(control$itmax[1]),
                     tol = as.double(control$tol[1]),
                     eps = as.double(control$eps),
                     niterout = integer(1), errout = double(1),
                     lwork = as.integer(lwork), info = FALSE)
        WARNING <- "EVV model is not available with prior"
        if(warn) warning(WARNING)
        temp <- structure(temp, info = NA, WARNING = WARNING, returnCode = -1)
        return(temp)
    }
    #
    z <- matrix(temp$z, n,K)
    loglik <- temp$loglik
    mu <- matrix(temp$mu, p,G)
    dimnames(mu) <- list(NULL, as.character(1:G))
    scale <- temp$scale[1]
    shape <- matrix(temp$shape, p,G)
    O <- aperm( array(temp$O, c(p,p,G)), c(2,1,3) )
    shape.o <- array( temp$shape.o, c(p,p,G) )
    pro <- temp$pro
    niterout <- temp$niterout
    errout <- temp$errout
    lapackSVDinfo <- temp$info
    WARNING <- NULL
    if(!is.finite(loglik) | any(is.nan(scale)) |
       any(is.nan(shape)) | any(is.nan(O)))
      { loglik <- .Machine$double.xmax }
    
    #
    if(lapackSVDinfo) {
        if(lapackSVDinfo > 0) {
            WARNING <- "LAPACK DGESVD fails to converge"
        }
        else {
            WARNING <- "input error for LAPACK DGESVD"
        }
        if(warn) warning(WARNING)
        z[] <- O[] <- shape[] <- NA
        scale <- loglik <- NA
        sigma <- array(NA, c(p, p, G))
        ret <- -9
        #
    } else if( loglik > signif(.Machine$double.xmax, 6)) {
        WARNING <- "singular covariance"
        if(warn) warning(WARNING)
        shape[] <- NA
        mu[] <- pro[] <- z[] <- loglik <- NA
        sigma <- array(NA, c(p, p, G))
        ret <- -1
    } else if(loglik <  - signif(.Machine$double.xmax, 6)) {
        if(control$equalPro) {
            WARNING <- "a z column sum fell below threshold"
            if(warn) warning(WARNING)
        } else {
            WARNING <- "mixing proportion fell below threshold"
            if(warn) warning(WARNING)
        }
        mu[] <- pro[] <- z[] <- loglik <- NA
        sigma <- array(NA, c(p, p, G))
        ret <- if(control$equalPro) -2 else -3
    } else {
        # scale <- sum(scale)/n
        sigma <- scale * shape.o
        if(niterout >= control$itmax[1]) {
            WARNING <- "iteration limit reached"
            if(warn) warning(WARNING)
            niterout <-  - niterout
            ret <- 1
        }
        else ret <- 0
    }
    info <- list(iterations = niterout, error = errout)
    # info <- c(iterations = its, error = err)
    dimnames(z) <- list(dimnames(data)[[1]], NULL)
    dimnames(mu) <- list(dimnames(data)[[2]], NULL)
    dimnames(O) <- list(dimnames(data)[[2]], dimnames(data)[[2]], NULL)
    ## Sigma = scale * O %*% diag(shape) %*% t(O)
    variance <- list(modelName = "EVV", d = p, G = G, sigma = sigma,
                     scale = scale, shape = shape, orientation = O)
    parameters <- list(pro=pro, mean=mu, variance=variance, Vinv=Vinv)
    structure(list(modelName = "EVV", prior = prior, n = n, d = p, G = G,
                   z = z, parameters = parameters, control = control,
                   loglik = loglik),
              info = info, WARNING = WARNING, returnCode = ret)
}

####
mstepEVV <- function(data, z, prior = NULL, warn = NULL, ...)
{
    if(is.null(warn)) warn <- mclust.options("warn")
    dimdat <- dim(data)
    oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
    if(oneD || length(dimdat) != 2)
        stop("data should be a matrix or a vector")
    data <- as.matrix(data)
    n <- nrow(data)
    p <- ncol(data)
    ##
    z <- as.matrix(z)
    dimz <- dim(z)
    if(dimz[1] != n)
        stop("row dimension of z should equal data length")
    G <- dimz[2]
    if(all(is.na(z))) {
        WARNING <- "z is missing"
        if(warn) warning(WARNING)
        variance <- list(modelName = "EVV", d = p, G = G,
                         scale = NA, shape = rep(NA,p), orientation=array(NA,c(p,p,G)))
        parameters <- list(pro=rep(NA,G), mean=matrix(NA,p,G),
                           variance=variance)
        return(structure(list(modelName="EVV", prior=prior, n=n, d=p,
                              G=G, z=z, parameters=parameters),
                         WARNING = WARNING, returnCode = 9))
        
    }
    if(any(is.na(z)) || any(z < 0) || any(z > 1)) stop("improper specification of z")
    lwork <- max(3 * min(n, p) + max(n, p), 5 * min(n, p), G)
    #
    # MICHAEL from here-------------------------------------------------------
    #
    # without prior specification
    if(is.null(prior)) {
        temp <- .Fortran( "msevv",
                          x = as.double(data),
                          z = as.double(z),
                          n = as.integer(n),
                          p = as.integer(p),
                          G = as.integer(G),
                          mu = double(p*G),
                          O = double(p*p*G),
                          shape.o = double(p*p*G),
                          scale = double(G),
                          shape = double(p*G),
                          pro = double(G),
                          lwork = as.integer(lwork),
                          info = as.integer(0),
                          eps = as.double(.Machine$double.eps),
                          package = "mclust")
    } else {
        # with prior
        priorParams <- do.call(prior$functionName,
                               c(list(data = data, G = G, modelName = "EVV"),
                                 prior[names(prior) != "functionName"]))
        # temp <- .Fortran("meevvp", ...)
        temp <- list(x = data, z = z, n = n, p = p, G = G,
                     mu = double(p*G), O = double(p*p*G), shape.o = double(p*p*G),
                     scale = double(G), shape = double(p*G), pro = double(G),
                     lwork = as.integer(lwork), info = FALSE,
                     eps = as.double(.Machine$double.eps))
        WARNING <- "EVV model is not available with prior"
        if(warn) warning(WARNING)
    }
    #
    lapackSVDinfo <- temp$info
    mu <- matrix(temp$mu, p,G)
    dimnames(mu) <- list(NULL, as.character(1:G))
    scale <- temp$scale[1]                                  # lambda
    O <- aperm( array(temp$O, c(p,p,G)), c(2,1,3) )
    shape.o <- array( temp$shape.o, c(p,p,G) )
    shape <- matrix(temp$shape, p,G)
    pro <- temp$pro
    WARNING <- NULL
    #
    if(lapackSVDinfo) {
        if(lapackSVDinfo > 0) {
            WARNING <- "LAPACK DGESVD fails to converge"
            ret <- -4
        }
        else {
            WARNING <- "input error for LAPACK DGESVD"
            ret <- -5
        }
        if(warn) warning(WARNING)
        O[] <- shape[] <- scale[] <- NA
        sigma <- array(NA, c(p, p, G))
        #
    } else if( any(abs(c(scale, shape)) > signif(.Machine$double.xmax, 6)) ) {
        WARNING <- "cannot compute M-step"
        if(warn) warning(WARNING)
        mu[] <- pro[] <- scale <- shape[] <- O[] <- NA
        sigma <- array(NA, c(p, p, G))
        ret <- -1
    } else {
        #        scale <- sum(scale)/n      
        #        scale <- sum(scale)/sum(z)       # lambda --> if noise, see help(mstep)
        sigma <- scale * shape.o
        ret <- 0
    }
    dimnames(z) <- list(dimnames(data)[[1]], NULL)
    dimnames(mu) <- dimnames(shape) <- list(dimnames(data)[[2]], NULL)
    dimnames(sigma) <- dimnames(O) <-
        list(dimnames(data)[[2]], dimnames(data)[[2]], NULL)
    variance <- list(modelName = "EVV", d = p, G = G, sigma = sigma,
                     scale = scale, shape = shape, orientation = O)
    parameters <- list(pro=pro, mean=mu, variance=variance)
    structure(list(modelName = "EVV", prior = prior, n = n, d = p, G = G,
                   z = z, parameters = parameters),
              WARNING = WARNING, returnCode = ret)
}

####
estepEVV <- function(data, parameters, warn = NULL, ...)
{
    if(is.null(warn)) warn <- mclust.options("warn")
    dimdat <- dim(data)
    if(is.null(dimdat) || length(dimdat) != 2)
        stop("data must be a matrix")
    data <- as.matrix(data)
    n <- nrow(data)
    p <- ncol(data)
    pro <- parameters$pro
    pro <- pro/sum(pro)
    l <- length(pro)
    mu <- as.matrix(parameters$mean)
    scale <- parameters$variance$scale
    shape <- parameters$variance$shape
    O <- parameters$variance$orientation
    G <- ncol(mu)
    noise <- l == G + 1
    if(!noise) {
        if(l != G)
            stop("pro improperly specified")
        K <- G
        Vinv <- NULL
    } else {
        K <- G + 1
        Vinv <- parameters$Vinv
        if(is.null(Vinv) || Vinv <= 0)
            Vinv <- hypvol(data, reciprocal = TRUE)
    }
    if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
           any(is.null(parameters[c("pro", "mean", "variance")]))) {
        WARNING <- "parameters are missing"
        if(warn) warning(WARNING)
        z <- matrix(NA,n,K)
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        return(structure(list(modelName = "EVV", n=n, d=p, G=G, z=z,
                              parameters=parameters, loglik=NA),
                         WARNING = WARNING, returnCode = 9))
    }
    if (is.null(parameters$variance$scale) ||
            is.null(parameters$variance$shape) ||
            is.null(parameters$variance$orientation))
        stop("variance parameters are missing")
    #
    # MICHAEL from here----------------------------------------------
    #
    temp <- .Fortran( "esevv",
                      x = as.double(data),
                      z = double(n*K),
                      n = as.integer(n),
                      p = as.integer(p),
                      G = as.integer(G),
                      Gnoise = as.integer(K),
                      mu = as.double(mu),
                      O =  as.double( aperm(O, c(2,1,3)) ),
                      scale = as.double(scale),
                      shape = as.double(shape),
                      pro = as.double(pro),
                      Vinv = as.double( if (is.null(Vinv)) -1 else Vinv ),
                      loglik = double(1),
                      eps = as.double(.Machine$double.eps),
                      package = "mclust")
    #
    loglik <- temp$loglik
    z <- matrix(temp$z, n,K)
    WARNING <- NULL
    if(loglik > signif(.Machine$double.xmax, 6)) {
        WARNING <- "singular covariance"
        if(warn) warning(WARNING)
        z[] <- loglik <- NA
        ret <- -1
    }
    else ret <- 0
    dimnames(z) <- list(dimnames(data)[[1]],NULL)
    structure(list(modelName = "EVV", n = n, d = p, G = G,
                   z = z, parameters = parameters, loglik = loglik),
              WARNING = WARNING, returnCode = ret)
}

####
cdensEVV <- function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
    if(is.null(warn)) warn <- mclust.options("warn")
    dimdat <- dim(data)
    if(is.null(dimdat) || length(dimdat) != 2)
        stop("data must be a matrix")
    data <- as.matrix(data)
    n <- nrow(data)
    p <- ncol(data)
    pro <- parameters$pro
    pro <- pro/sum(pro)
    l <- length(pro)
    mu <- as.matrix(parameters$mean)
    scale <- parameters$variance$scale
    shape <- parameters$variance$shape
    O <- parameters$variance$orientation
    p <- ncol(data)
    G <- ncol(mu)
    if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
           any(is.null(parameters[c("pro", "mean", "variance")]))) {
        WARNING <- "parameters are missing"
        if(warn) warning(WARNING)
        z <- matrix(NA,n,G)
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        return(structure(z, logarithm = logarithm, modelName = "EVV",
                         WARNING = WARNING, returnCode = 9))
    }
    if (is.null(parameters$variance$scale) ||
            is.null(parameters$variance$shape) ||
            is.null(parameters$variance$orientation))
        stop("variance parameters are missing")
    #
    temp <- .Fortran( "esevv",
                      x = as.double(data),
                      z = double(n*G),
                      n = as.integer(n),
                      p = as.integer(p),
                      G = as.integer(G),
		      Gnoise = as.integer(G),
                      mu = as.double(mu),
                      O =  as.double( aperm(O, c(2,1,3)) ),
                      scale = as.double(scale),
                      shape = as.double(shape),
                      pro = as.double(-1),
		      Vinv = as.double(-1),
                      loglik = double(1),
                      eps = as.double(.Machine$double.eps),
                      package = "mclust")
    #
    loglik <- temp$loglik
    z <- matrix(temp$z, n,G)
    WARNING <- NULL
    if(loglik > signif(.Machine$double.xmax, 6)) {
        WARNING <- "singular covariance"
        if(warn) warning(WARNING)
        z[] <- NA
        ret <- -1
    } else {
        if (!logarithm) z <- exp(z)
        ret <- 0
    }
    dimnames(z) <- list(dimnames(data)[[1]],NULL)
    structure(z, logarithm = logarithm, modelName = "EVV",
              WARNING = WARNING, returnCode = ret)
}

###
simEVV <- function(parameters, n, seed = NULL, ...) 
{
    if (!is.null(seed)) 
        set.seed(seed)
    mu <- as.matrix(parameters$mean)
    d <- nrow(mu)
    G <- ncol(mu)
    if (any(is.na(parameters[c("mean", "variance")])) || 
            any(is.null(parameters[c("mean", "variance")]))) {
        warn <- "parameters are missing"
        warning("parameters are missing")
        return(structure(matrix(as.double(NA), n, d + 1), modelName = "EVV"))
    }
    pro <- parameters$pro
    if (is.null(pro)) 
        pro <- rep(1/G, G)
    clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
    ctabel <- tabulate(clabels, nbins = G)
    x <- matrix(0, n, d)
    rtshape <- sqrt(parameters$variance$shape)
    if (dim(rtshape)[1] != d | dim(rtshape)[2] != G) 
        stop("shape incompatible with mean")
    rtscale <- sqrt(parameters$variance$scale)    
    for (k in 1:G) 
    {
        m <- ctabel[k]
        sss <- rtscale * rtshape[,k]
        cholSigma <- t(parameters$variance$orientation[,,k]) * sss
        x[clabels == k, ] <- sweep( matrix(rnorm(m*d), nrow = m, ncol = d) %*% cholSigma, 
                                   MARGIN = 2, STATS = mu[,k], FUN = "+" )
    }
    dimnames(x) <- list(NULL, 1:d)
    structure(cbind(group = clabels, x), modelName = "EVV")
}


##############################################################################
###                               VEE model                               ####
##############################################################################
emVEE <- function(data, parameters, prior = NULL, control = emControl(),
                  warn = NULL, ...)
{
    z <- estepVEE(data, parameters = parameters, warn = warn)$z
    meVEE(data, z = z, prior = prior, control = control,
          Vinv = parameters$Vinv, warn = warn)
}

####
meVEE <- function(data, z, prior = NULL, control = emControl(),
                  Vinv = NULL, warn = NULL, ...)
{
    if(is.null(warn)) warn <- mclust.options("warn")
    dimdat <- dim(data)
    oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
    if(oneD || length(dimdat) != 2)
        stop("data should in the form of a matrix")
    data <- as.matrix(data)
    n <- nrow(data)
    p <- ncol(data)
    z <- as.matrix(z)
    dimz <- dim(z)
    if(dimz[1] != n)
        stop("data and z should have the same row dimension")
    K <- dimz[2]
    if (!is.null(Vinv)) {
        G <- K - 1
        if(Vinv <= 0) Vinv <- hypvol(data, reciprocal = TRUE)
    } else G <- K
    if(all(is.na(z))) {
        WARNING <- "z is missing"
        if(warn) warning(WARNING)
        variance <- list(modelName = "VEE", d = p, G = G,
                         scale=rep(NA,G), shape=rep(NA,p), orientation=array(NA,c(p,p,G)))
        parameters <- list(pro=rep(NA,G), mean=matrix(NA,p,G),
                           variance=variance)
        return(structure(list(modelName="VEE", prior=prior, n=n, d=p,
                              G=G, z=z, parameters=parameters,
                              control=control, loglik=NA),
                         WARNING = WARNING, returnCode = 9))
    }
    if(any(is.na(z)) || any(z < 0) || any(z > 1))
        stop("improper specification of z")
    lwork <- max(3 * min(n, p) + max(n, p), 5 * min(n, p), p + G)
    storage.mode(z) <- "double"
    #
    # MICHAEL from here-------------------------------------------------------
    #
    # without prior specification
    if(is.null(prior)) {
        temp <- .Fortran("mevee",
                         x = as.double(data),
                         z = as.double(z),
                         n = as.integer(n),
                         p = as.integer(p),
                         G = as.integer(G),
                         Gnoise = as.integer(K),
                         mu = double(p*G),
                         C = double(p*p),
                         U = double(p*p*G),
                         scale = double(G),
                         shape = double(p),
                         pro = double(K),
                         Vinv = as.double( if (is.null(Vinv)) -1 else Vinv ),
                         loglik = double(1),
                         eqpro = as.logical(control$equalPro),
                         itmaxin = as.integer(control$itmax[2]),
                         tolin = as.double(control$tol[2]),
                         itmaxout = as.integer(control$itmax[1]),
                         tolout = as.double(control$tol[1]),
                         eps = as.double(control$eps),
                         niterin = integer(1),
                         errin = double(1),
                         niterout = integer(1),
                         errout = double(1),
                         lwork = as.integer(lwork),
                         info = as.integer(0),
                         package = "mclust")
        #
    } else {
        # with prior
        priorParams <- do.call(prior$functionName,
                               c(list(data = data, G = G, modelName = "VEE"),
                                 prior[names(prior) != "functionName"]))
        # temp <- .Fortran("meveep", ...)
        temp <- list(x = data, z = z, n = n, p = p, G = G,
                     mu = double(p*G), C = double(p*p), U = double(p*p*G),
                     scale = double(G), shape = double(p), pro = double(K),
                     loglik = NA, eqpro = as.logical(control$equalPro),
                     itmaxin = as.integer(control$itmax[2]),
                     tolin = as.double(control$tol[2]),
                     itmaxout = as.integer(control$itmax[1]),
                     tolout = as.double(control$tol[1]),
                     eps = as.double(control$eps),
                     niterin = integer(1), errin = double(1),
                     niterout = integer(1), errout = double(1),
                     lwork = as.integer(lwork), info = FALSE)
        WARNING <- "VEE model is not available with prior"
        if(warn) warning(WARNING)
        temp <- structure(temp, info = NA, WARNING = WARNING, returnCode = -1)
        return(temp)
    }
    z <- matrix(temp$z, n,K)
    niterin <- temp$niterin
    errin <- temp$errin
    niterout <- temp$niterout
    errout <- temp$errout
    loglik <- temp$loglik
    lapackSVDinfo <- temp$info
    mu <- matrix(temp$mu, p,G)
    dimnames(mu) <- list(NULL, as.character(1:G))
    scale <- temp$scale
    shape <- temp$shape
    shape.o <- matrix(temp$C, p,p)
    O <- if(any(is.nan(shape.o))) shape.o else
         svd(shape.o, nu = 0)$v
    pro <- temp$pro
    WARNING <- NULL
    if(!is.finite(loglik) | any(is.nan(scale)) |
       any(is.nan(shape)) | any(is.nan(O)))
      { loglik <- .Machine$double.xmax }
    #
    if(lapackSVDinfo) {
        if(lapackSVDinfo > 0) {
            WARNING <- "LAPACK DSYEV or DPOTRI fails to converge"
        }
        else {
            WARNING <- "input error for LAPACK DPOTRF, DSYEV or DPOTRI"
        }
        if(warn) warning(WARNING)
        O[] <- shape[] <- scale[] <- NA
        mu[] <- pro[] <- z[] <- loglik <- NA
        sigma <- array(NA, c(p, p, G))
        ret <- -9
        #
    } else if(loglik > signif(.Machine$double.xmax, 6)) {
        WARNING <- "singular covariance"
        if(warn) warning(WARNING)
        O[] <- shape[] <- scale[] <- NA
        mu[] <- pro[] <- z[] <- loglik <- NA
        sigma <- array(NA, c(p, p, G))
        ret <- -1
    } else if(loglik <  - signif(.Machine$double.xmax, 6)) {
        if(control$equalPro) {
            WARNING <- "z column sum fell below threshold"
            if(warn) warning(WARNING)
        } else {
            WARNING <- "mixing proportion fell below threshold"
            if(warn) warning(WARNING)
        }
        mu[] <- pro[] <- z[] <- loglik <- NA
        sigma <- array(NA, c(p, p, G))
        ret <- if(control$equalPro) -2 else -3
    } else {
        sigma <- sweep( array(shape.o, c(p,p,G)), 3, FUN = "*", STATS = scale )
        if(niterin >= control$itmax[2]) {
            WARNING <- "inner iteration limit reached"
            if(warn) warning(WARNING)
            niterin <-  - niterin
            ret <- 2
        } else if(niterout >= control$itmax[1]) {
            WARNING <- "iteration limit reached"
            if(warn) warning(WARNING)
            niterout <-  - niterout
            ret <- 1
        } else ret <- 0
    }
    info <- structure(c(niterout = niterout, errout = errout),
                      inner = c(niterin = niterin, errin = errin))
    # info <- structure(c(iterations = its, error = err),
    #                   inner = c(iterations = inner, error = inerr))
    dimnames(z) <- list(dimnames(data)[[1]],NULL)
    dimnames(mu) <- list(dimnames(data)[[2]], NULL)
    dimnames(sigma) <- dimnames(O) <- list(dimnames(data)[[2]], dimnames(data)[[2]])
    ##  Sigma = scale * O %*% diag(shape) %*% t(O)
    variance <- list(modelName = "VEE", d = p, G = G, sigma = sigma,
                     scale = scale, shape = shape, orientation = O)
    parameters <- list(pro=pro, mean=mu, variance=variance, Vinv=Vinv)
    structure(list(modelName = "VEE", prior = prior, n = n, d = p, G = G,
                   z = z, parameters = parameters, control = control,
                   loglik = loglik),
              info = info, WARNING = WARNING, returnCode = ret)
}

####
mstepVEE <- function(data, z, prior = NULL, warn = NULL, control = NULL, ...)
{
    if(is.null(warn)) warn <- mclust.options("warn")
    dimdat <- dim(data)
    oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
    if(oneD || length(dimdat) != 2)
        stop("data should be a matrix or a vector")
    data <- as.matrix(data)
    n <- nrow(data)
    p <- ncol(data)
    ##
    z <- as.matrix(z)
    dimz <- dim(z)
    if(dimz[1] != n)
        stop("row dimension of z should equal data length")
    G <- dimz[2]
    if(all(is.na(z))) {
        WARNING <- "z is missing"
        if(warn) warning(WARNING)
        variance <- list(modelName = "VEE", d = p, G = G,
                         scale = rep(NA,G), shape = rep(NA,p), orientation = array(NA,c(p,p,G)))
        parameters <- list(pro=rep(NA,G), mean=matrix(NA,p,G),
                           variance=variance)
        return(structure(list(modelName="VEE", prior=prior, n=n, d=p,
                              G=G, z=z, parameters=parameters,
                              control=control, loglik=NA),
                         WARNING = WARNING, returnCode = 9))
        
        WARNING <- "z is missing"
        if(warn) warning(WARNING)
        return(structure(list(
            n = n, d = p, G = G, mu = matrix(NA,p, G), sigma = array(NA, c(p, p, G)),
            decomp = list(d = p, G = G, scale = rep(NA, G), shape = rep(NA, p),
                          orientation = array(NA, c(p, p, G))),
            pro = rep(NA,G), modelName = "VEE", prior = prior), WARNING = WARNING))
    }
    #  shape <- sqrt(rev(sort(shape/exp(sum(log(shape))/p))))
    if(any(is.na(z)) || any(z < 0) || any(z > 1)) stop("improper specification of z")
    if (is.null(control)) control <- emControl()
    itmax <- if(length(control$itmax) == 1) control$itmax else control$itmax[2]
    tol <- if(length(control$tol) == 1) control$tol else control$tol[2]
    lwork <- max(3 * min(n, p) + max(n, p), 5 * min(n, p), p + G)
    #
    # MICHAEL from here-------------------------------------------------------
    #
    # without prior specification
    if(is.null(prior)) {
        temp <- .Fortran( "msvee",
                          x = as.double(data),
                          z = as.double(z),
                          n = as.integer(n),
                          p = as.integer(p),
                          G = as.integer(G),
                          mu = double(p*G),
                          U = double(p*p*G),
                          C = double(p*p),
                          scale = as.double( rep(1,G) ),
                          pro = double(G),
                          lwork = as.integer(lwork),
                          info = as.integer(0),
                          itmax = as.integer(itmax),
                          tol = as.double(tol),
                          niterin = integer(1),
                          errin = double(1),
                          eps = as.double(.Machine$double.eps),
                          package = "mclust")
    } else {
        # with prior
        priorParams <- do.call(prior$functionName,
                               c(list(data = data, G = G, modelName = "VEE"),
                                 prior[names(prior) != "functionName"]))
        #
        # temp <- .Fortran("msveep", ...)
        temp <- list(x = data, z = z, n = n, p = p, G = G,
                     mu = double(p*G),  U = double(p*p*G), C = double(p*p),
                     scale = double(G), pro = double(G),
                     lwork = as.integer(lwork), info = FALSE,
                     itmax = as.integer(itmax), tol = as.double(tol),
                     niterin = integer(1), errin = double(1),
                     eps = as.double(control$eps))
        WARNING <- "VEE model is not available with prior"
        if(warn) warning(WARNING)
    }
    lapackSVDinfo <- temp$info
    errin <- temp$errin
    niterin <- temp$niterin
    mu <- matrix(temp$mu, p,G)
    dimnames(mu) <- list(NULL, as.character(1:G))
    scale <- temp$scale
    shape.o <- matrix(temp$C, p,p)
    SVD <- svd(shape.o, nu = 0)
    shape <- SVD$d
    O <- SVD$v
    pro <- temp$pro
    WARNING <- NULL
    #
    if(lapackSVDinfo) {
        if(lapackSVDinfo > 0) {
            WARNING <- "LAPACK DSYEV or DPOTRI fails to converge"
        } else {
            WARNING <- "input error for LAPACK DPOTRF, DSYEV or DPOTRI"
        }
        if(warn) warning(WARNING)
        O[] <- shape[] <- scale[] <- NA
        sigma <- array(NA, c(p, p, G))
        ret <- -9
        #
    } else if(any(c(scale, shape) > signif(.Machine$double.xmax, 6))) {
        WARNING <- "cannot compute M-step"
        if(warn) warning(WARNING)
        mu[] <- pro[] <- O[] <- shape[] <- scale[] <- NA
        sigma <- array(NA, c(p, p, G))
        ret <- -1
    } else {
        sigma <- sweep( array(shape.o, c(p,p,G)), 3, FUN = "*", STATS = scale )
        if(niterin >= itmax) {
            WARNING <- "inner iteration limit reached"
            if(warn) warning(WARNING)
            niterin <-  - niterin
        }
        ret <- 2
    }
    info <- c(iteration = niterin, error = errin)
    dimnames(z) <- list(dimnames(data)[[1]], NULL)
    dimnames(mu) <- list(dimnames(data)[[2]], NULL)
    dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]])
    dimnames(O) <- list(dimnames(data)[[2]], dimnames(data)[[2]])
    variance <- list(modelName = "VEE", d = p, G = G, sigma = sigma,
                     scale = scale, shape = shape, orientation = O)
    parameters <- list(pro=pro, mean=mu, variance=variance)
    structure(list(modelName = "VEE", prior = prior, n = n, d = p, G = G,
                   z = z, parameters = parameters, control = control),
              info = info, WARNING = WARNING, returnCode = ret)
}

###
estepVEE <- function(data, parameters, warn = NULL, ...)
{
    if(is.null(warn)) warn <- mclust.options("warn")
    dimdat <- dim(data)
    if(is.null(dimdat) || length(dimdat) != 2)
        stop("data must be a matrix")
    data <- as.matrix(data)
    n <- nrow(data)
    p <- ncol(data)
    pro <- parameters$pro
    pro <- pro/sum(pro)
    l <- length(pro)
    mu <- as.matrix(parameters$mean)
    scale <- parameters$variance$scale
    shape <- parameters$variance$shape
    O <- parameters$variance$orientation
    G <- ncol(mu)
    noise <- l == G + 1
    if(!noise) {
        if(l != G)
            stop("pro improperly specified")
        K <- G
        Vinv <- NULL
    } else {
        K <- G + 1
        if(is.null(Vinv) || Vinv <= 0)
            Vinv <- hypvol(data, reciprocal = TRUE)
    }
    if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
           any(is.null(parameters[c("pro", "mean", "variance")]))) {
        WARNING <- "parameters are missing"
        if(warn) warning(WARNING)
        z <- matrix(NA,n,K)
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        return(structure(list(modelName = "VEE", n=n, d=p, G=G, z=z,
                              parameters=parameters, loglik=NA),
                         WARNING = WARNING, returnCode = 9))
    }
    if (is.null(parameters$variance$scale) ||
            is.null(parameters$variance$shape) ||
            is.null(parameters$variance$orientation))
        stop("variance parameters are missing")
    #
    # MICHAEL from here-------------------------------------------------------
    #
    temp <- .Fortran( "esvee",
                      x = as.double(data),
                      z = double(n*K),
                      n = as.integer(n),
                      p = as.integer(p),
                      G = as.integer(G),
                      Gnoise = as.integer(K),
                      mu = as.double(mu),
                      O =  as.double( t(O) ),
                      scale = as.double(scale),
                      shape = as.double(shape),
                      pro = as.double(pro),
		      Vinv = as.double( if (is.null(Vinv)) -1 else Vinv ),
                      loglik = double(1),
                      eps = as.double(.Machine$double.eps),
                      package = "mclust")
    #
    loglik <- temp$loglik
    z <- matrix(temp$z, n,K)
    WARNING <- NULL
    if(loglik > signif(.Machine$double.xmax, 6)) {
        WARNING <- "cannot compute E-step"
        if(warn) warning(WARNING)
        z[] <- loglik <- NA
        ret <- -1
    }
    else ret <- 0
    dimnames(z) <- list(dimnames(data)[[1]],NULL)
    structure(list(modelName = "VEE", n = n, d = p, G = G,
                   z = z, parameters = parameters, loglik = loglik),
              WARNING = WARNING, returnCode = ret)
    
}

####
cdensVEE <- function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
    if(is.null(warn)) warn <- mclust.options("warn")
    dimdat <- dim(data)
    if(is.null(dimdat) || length(dimdat) != 2)
        stop("data must be a matrix")
    data <- as.matrix(data)
    n <- nrow(data)
    p <- ncol(data)
    pro <- parameters$pro
    pro <- pro/sum(pro)
    l <- length(pro)
    mu <- as.matrix(parameters$mean)
    scale <- parameters$variance$scale
    shape <- parameters$variance$shape
    O <- parameters$variance$orientation
    p <- ncol(data)
    G <- ncol(mu)
    if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
           any(is.null(parameters[c("pro", "mean", "variance")]))) {
        WARNING <- "parameters are missing"
        if(warn) warning(WARNING)
        z <- matrix(NA,n,G)
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        return(structure(z, logarithm = logarithm, modelName = "VEE",
                         WARNING = WARNING, returnCode = 9))
    }
    if (is.null(parameters$variance$scale) ||
            is.null(parameters$variance$shape) ||
            is.null(parameters$variance$orientation))
        stop("variance parameters are missing")
    #
    # MICHAEL from here-------------------------------------------------------
    #
    temp <- .Fortran( "esvee",
                      x = as.double(data),
                      z = double(n*G),
                      n = as.integer(n),
                      p = as.integer(p),
                      G = as.integer(G),
		      Gnoise = as.integer(G),
                      mu = as.double(mu),
                      O =  as.double( t(O) ),
                      scale = as.double(scale),
                      shape = as.double(shape),
                      pro = as.double(-1),
		      Vinv = as.double(-1),
                      loglik = double(1),
                      eps = as.double(.Machine$double.eps),
                      package = "mclust")
    #
    loglik <- temp$loglik
    z <- matrix(temp$z, n,G)
    WARNING <- NULL
    if(loglik > signif(.Machine$double.xmax, 6)) {
        WARNING <- "cannot compute E-step"
        if(warn) warning(WARNING)
        z[] <- NA
        ret <- -1
    } else {
        if (!logarithm) z <- exp(z)
        ret <- 0
    }
    dimnames(z) <- list(dimnames(data)[[1]],NULL)
    structure(z, logarithm = logarithm, modelName = "VEE",
              WARNING = WARNING, returnCode = ret)
}

###
simVEE <- function(parameters, n, seed = NULL, ...) 
{
    if (!is.null(seed)) 
        set.seed(seed)
    mu <- as.matrix(parameters$mean)
    d <- nrow(mu)
    G <- ncol(mu)
    if (any(is.na(parameters[c("mean", "variance")])) || 
            any(is.null(parameters[c("mean", "variance")]))) {
        warn <- "parameters are missing"
        warning("parameters are missing")
        return(structure(matrix(as.double(NA), n, d + 1), modelName = "VEE"))
    }
    pro <- parameters$pro
    if (is.null(pro)) 
        pro <- rep(1/G, G)
    clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
    ctabel <- tabulate(clabels, nbins = G)
    x <- matrix(0, n, d)
    rtshape <- sqrt(parameters$variance$shape)
    if (length(rtshape) != d) 
        stop("shape incompatible with mean")
    rtscale <- sqrt(parameters$variance$scale)    
    if (length(rtscale) != G) 
        stop("scale incompatible with mean")
    for (k in 1:G) 
    {
        m <- ctabel[k]
        sss <- rtscale[k] * rtshape
        cholSigma <- t(parameters$variance$orientation) * sss
        x[clabels == k, ] <- sweep( matrix(rnorm(m*d), nrow = m, ncol = d) %*% cholSigma, 
                                    MARGIN = 2, STATS = mu[,k], FUN = "+" )
    }
    dimnames(x) <- list(NULL, 1:d)
    structure(cbind(group = clabels, x), modelName = "VEE")
}


##############################################################################
###                               EVE model                               ####
##############################################################################
emEVE <- function(data, parameters, prior = NULL, control = emControl(),
                  warn = NULL, ...)
{
    z <- estepEVE(data, parameters = parameters, warn = warn)$z
    meEVE(data, z = z, prior = prior, control = control,
          Vinv = parameters$Vinv, warn = warn)
}

####
meEVE <- function(data, z, prior = NULL, control = emControl(),
                  Vinv = NULL, warn = NULL, ...)
{
    if(is.null(warn)) warn <- mclust.options("warn")
    dimdat <- dim(data)
    oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
    if(oneD || length(dimdat) != 2)
        stop("data should in the form of a matrix")
    data <- as.matrix(data)
    n <- nrow(data)
    p <- ncol(data)
    z <- as.matrix(z)
    dimz <- dim(z)
    if(dimz[1] != n)
        stop("data and z should have the same row dimension")
    K <- dimz[2]
    if (!is.null(Vinv)) {
        G <- K - 1
        if(Vinv <= 0) Vinv <- hypvol(data, reciprocal = TRUE)
    } else G <- K
    if(all(is.na(z))) {
        WARNING <- "z is missing"
        if(warn) warning(WARNING)
        variance <- list(modelName = "EVE", d = p, G = G,
                         scale=rep(NA,G), shape=rep(NA,p), orientation=array(NA,c(p,p,G)))
        parameters <- list(pro=rep(NA,G), mean=matrix(NA,p,G),
                           variance=variance)
        return(structure(list(modelName="EVE", prior=prior, n=n, d=p,
                              G=G, z=z, parameters=parameters,
                              control=control, loglik=NA),
                         WARNING = WARNING, returnCode = 9))
    }
    if(any(is.na(z)) || any(z < 0) || any(z > 1))
        stop("improper specification of z")
    lwork <- max(3 * min(n, p) + max(n, p), 5 * min(n, p), p + G)
    storage.mode(z) <- "double"
    
    #
    # MICHAEL from here-------------------------------------------------------
    #
    # without prior specification
    if(is.null(prior)) {
        temp <- .Fortran("meeve",
                         x = as.double(data),
                         z = as.double(z),
                         n = as.integer(n),
                         p = as.integer(p),
                         G = as.integer(G),
                         Gnoise = as.integer(K),
                         mu = double(p*G),
                         O = as.double( diag(p) ),
                         U = double(p*p*G),
                         scale = double(1),
                         shape = as.double( matrix(1, p,G) ),
                         pro = double(K),
                         Vinv = as.double( if (is.null(Vinv)) -1 else Vinv ),
                         loglik = double(1),
                         eqpro = as.logical(control$equalPro),
                         itmaxin = as.integer(control$itmax[2]),
                         tolin = as.double(control$tol[2]),
                         itmaxout = as.integer(control$itmax[1]),
                         tolout = as.double(control$tol[1]),
                         eps = as.double(control$eps),
                         niterin = integer(1),
                         errin = double(1),
                         niterout = integer(1),
                         errout = double(1),
                         lwork = as.integer(lwork),
                         info = as.integer(0),
                         package = "mclust")
        #
    } else {
        # with prior
        priorParams <- do.call(prior$functionName,
                               c(list(data = data, G = G, modelName = "EVE"),
                                 prior[names(prior) != "functionName"]))
        # temp <- .Fortran("meevep", ...)
        temp <- list(x = data, z = z, n = n, p = p, G = G,
                     mu = double(p*G), O = double(p*p), U = double(p*p*G),
                     scale = double(1), shape = double(p*G), pro = double(G),
                     loglik = NA, eqpro = as.logical(control$equalPro),
                     itmaxin = as.integer(control$itmax[2]),
                     tolin = as.double(control$tol[2]),
                     itmaxout = as.integer(control$itmax[1]),
                     tolout = as.double(control$tol[1]),
                     eps = as.double(control$eps),
                     niterin = integer(1), errin = double(1),
                     niterout = integer(1), errout = double(1),
                     lwork = as.integer(lwork), info = FALSE)
        WARNING <- "EVE model is not available with prior"
        if(warn) warning(WARNING)
        temp <- structure(temp, info = NA, WARNING = WARNING, returnCode = -1)
        return(temp)
    }
    z <- matrix(temp$z, n,K)
    niterin <- temp$niterin
    errin <- temp$errin
    niterout <- temp$niterout
    errout <- temp$errout
    loglik <- temp$loglik
    lapackSVDinfo <- temp$info
    mu <- matrix(temp$mu, p,G)
    dimnames(mu) <- list(NULL, as.character(1:G))
    scale <- temp$scale
    shape <- matrix(temp$shape, p,G)
    O <- t( matrix(temp$O, p,p) )
    pro <- temp$pro
    WARNING <- NULL
    if(!is.finite(loglik) | any(is.nan(scale)) |
       any(is.nan(shape)) | any(is.nan(O)))
      { loglik <- .Machine$double.xmax }
    #
    if(lapackSVDinfo) {
        if(lapackSVDinfo > 0) {
            WARNING <- "LAPACK DSYEV or DGESVD fails to converge"
        }
        else {
            WARNING <- "input error for LAPACK DSYEV or DGESVD"
        }
        if(warn) warning(WARNING)
        O[] <- shape[] <- scale[] <- NA
        mu[] <- pro[] <- z[] <- loglik <- NA
        sigma <- array(NA, c(p, p, G))
        ret <- -9
        #
    } else if(loglik > signif(.Machine$double.xmax, 6)) {
        WARNING <- "singular covariance"
        if(warn) warning(WARNING)
        O[] <- shape[] <- scale[] <- NA
        mu[] <- pro[] <- z[] <- loglik <- NA
        sigma <- array(NA, c(p, p, G))
        ret <- -1
    } else if(loglik <  - signif(.Machine$double.xmax, 6)) {
        if(control$equalPro) {
            WARNING <- "z column sum fell below threshold"
            if(warn) warning(WARNING)
        } else {
            WARNING <- "mixing proportion fell below threshold"
            if(warn) warning(WARNING)
        }
        mu[] <- pro[] <- z[] <- loglik <- NA
        sigma <- array(NA, c(p, p, G))
        ret <- if(control$equalPro) -2 else -3
    } else {
        sigma <- array( apply(shape, 2, function(sh) scale * O%*%diag(sh)%*%t(O)), c(p,p,G) )
        if(niterin >= control$itmax[2]) {
            WARNING <- "inner iteration limit reached"
            if(warn) warning(WARNING)
            niterin <-  - niterin
            ret <- 2
        } else if(niterout >= control$itmax[1]) {
            WARNING <- "iteration limit reached"
            if(warn) warning(WARNING)
            niterout <-  - niterout
            ret <- 1
        } else ret <- 0
    }
    info <- structure(c(niterout = niterout, errout = errout),
                      inner = c(niterin = niterin, errin = errin))
    # info <- structure(c(iterations = its, error = err),
    #                   inner = c(iterations = inner, error = inerr))
    dimnames(z) <- list(dimnames(data)[[1]],NULL)
    dimnames(mu) <- list(dimnames(data)[[2]], NULL)
    dimnames(sigma) <- dimnames(O) <- list(dimnames(data)[[2]], dimnames(data)[[2]])
    ##  Sigma = scale * O %*% diag(shape) %*% t(O)
    variance <- list(modelName = "EVE", d = p, G = G, sigma = sigma,
                     scale = scale, shape = shape, orientation = O)
    parameters <- list(pro=pro, mean=mu, variance=variance, Vinv=Vinv)
    structure(list(modelName = "EVE", prior = prior, n = n, d = p, G = G,
                   z = z, parameters = parameters, control = control,
                   loglik = loglik),
              info = info, WARNING = WARNING, returnCode = ret)
}


####
mstepEVE <- function(data, z, prior = NULL, warn = NULL, control = NULL, ...)
{
    if(is.null(warn)) warn <- mclust.options("warn")
    dimdat <- dim(data)
    oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
    if(oneD || length(dimdat) != 2)
        stop("data should be a matrix or a vector")
    data <- as.matrix(data)
    n <- nrow(data)
    p <- ncol(data)
    ##
    z <- as.matrix(z)
    dimz <- dim(z)
    if(dimz[1] != n)
        stop("row dimension of z should equal data length")
    G <- dimz[2]
    if(all(is.na(z))) {
        WARNING <- "z is missing"
        if(warn) warning(WARNING)
        variance <- list(modelName = "EVE", d = p, G = G,
                         scale = rep(NA,G), shape = rep(NA,p), orientation = array(NA,c(p,p,G)))
        parameters <- list(pro=rep(NA,G), mean=matrix(NA,p,G),
                           variance=variance)
        return(structure(list(modelName="EVE", prior=prior, n=n, d=p,
                              G=G, z=z, parameters=parameters,
                              control=control, loglik=NA),
                         WARNING = WARNING, returnCode = 9))
        
        WARNING <- "z is missing"
        if(warn) warning(WARNING)
        return(structure(list(
            n = n, d = p, G = G, mu = matrix(NA,p, G), sigma = array(NA, c(p, p, G)),
            decomp = list(d = p, G = G, scale = rep(NA, G), shape = rep(NA, p),
                          orientation = array(NA, c(p, p, G))),
            pro = rep(NA,G), modelName = "EVE", prior = prior), WARNING = WARNING))
    }
    if(any(is.na(z)) || any(z < 0) || any(z > 1)) stop("improper specification of z")
    if (is.null(control)) control <- emControl()
    itmax <- if(length(control$itmax) == 1) control$itmax else control$itmax[2]
    tol <- if(length(control$tol) == 1) control$tol else control$tol[2]
    lwork <- max(3 * min(n, p) + max(n, p), 5 * min(n, p), p + G)
    
    #
    # MICHAEL from here-------------------------------------------------------
    #
    # without prior specification
    if(is.null(prior)) {
        temp <- .Fortran("mseve",
                         x = as.double(data), 
                         z = as.double(z),
                         n = as.integer(n),
                         p = as.integer(p),
                         G = as.integer(G),
                         mu = double(p*G),
                         U = double(p*p*G),
                         O = as.double( diag(p) ),
                         scale = as.double(1),
                         shape = as.double( matrix(1, p,G) ),
                         pro = double(G),
                         lwork = as.integer(lwork),
                         info = as.integer(0),
                         itmax = as.integer(itmax),
                         tol = as.double(tol),
                         niterin = integer(1),
                         errin = double(1),
                         eps = as.double(.Machine$double.eps),
                         # d = 100000,
                         # trgtvec = as.double(100000),
                         package = "mclust")
    } else {
        # with prior
        priorParams <- do.call(prior$functionName,
                               c(list(data = data, G = G, modelName = "EVE"),
                                 prior[names(prior) != "functionName"]))
        #
        # temp <- .Fortran("msevep", ...)
        temp <- list(x = data, z = z, n = n, p = p, G = G,
                     mu = double(p*G),  U = double(p*p*G), O = double(p*p),
                     scale = double(1), pro = double(G), shape = double(p*G),
                     lwork = as.integer(lwork), info = FALSE,
                     itmax = as.integer(itmax), tol = as.double(tol),
                     niterin = integer(1), errin = double(1),
                     eps = as.double(.Machine$double.eps))
        WARNING <- "EVE model is not available with prior"
        if(warn) warning(WARNING)
    }
    lapackSVDinfo <- temp$info
    errin <- temp$errin
    niterin <- temp$niterin
    mu <- matrix(temp$mu, p,G)
    dimnames(mu) <- list(NULL, as.character(1:G))
    scale <- temp$scale
    O <- t( matrix(temp$O, p,p) )
    shape <- matrix(temp$shape, p,G)
    pro <- temp$pro
    WARNING <- NULL
    #
    if(lapackSVDinfo) {
        if(lapackSVDinfo > 0) {
            WARNING <- "LAPACK DSYEV or DGESVD fails to converge"
        } else {
            WARNING <- "input error for LAPACK DSYEV or DGESVD"
        }
        if(warn) warning(WARNING)
        O[] <- shape[] <- scale[] <- NA
        sigma <- array(NA, c(p, p, G))
        ret <- -9
        #
    } else if( any(c(scale, shape) > signif(.Machine$double.xmax, 6)) ) {
        WARNING <- "cannot compute M-step"
        if(warn) warning(WARNING)
        mu[] <- pro[] <- O[] <- shape[] <- scale[] <- NA
        sigma <- array(NA, c(p, p, G))
        ret <- -1
    } else {
        sigma <- array( apply(shape, 2, function(sh) scale * O%*%diag(sh)%*%t(O)), c(p,p,G) )
        if(niterin >= itmax) {
            WARNING <- "inner iteration limit reached"
            if(warn) warning(WARNING)
            niterin <-  - niterin
        }
        ret <- 2
    }
    info <- c(iteration = niterin, error = errin)
    dimnames(z) <- list(dimnames(data)[[1]], NULL)
    dimnames(mu) <- list(dimnames(data)[[2]], NULL)
    dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]])
    dimnames(O) <- list(dimnames(data)[[2]], dimnames(data)[[2]])
    variance <- list(modelName = "EVE", d = p, G = G, sigma = sigma,
                     scale = scale, shape = shape, orientation = O)
    parameters <- list(pro=pro, mean=mu, variance=variance)
    structure(list(modelName = "EVE", prior = prior, n = n, d = p, G = G,
                   z = z, parameters = parameters, control = control),
              info = info, WARNING = WARNING, returnCode = ret)
}


###
estepEVE <- function(data, parameters, warn = NULL, ...)
{
    if(is.null(warn)) warn <- mclust.options("warn")
    dimdat <- dim(data)
    if(is.null(dimdat) || length(dimdat) != 2)
        stop("data must be a matrix")
    data <- as.matrix(data)
    n <- nrow(data)
    p <- ncol(data)
    pro <- parameters$pro
    pro <- pro/sum(pro)
    l <- length(pro)
    mu <- as.matrix(parameters$mean)
    scale <- parameters$variance$scale
    shape <- parameters$variance$shape
    O <- parameters$variance$orientation
    G <- ncol(mu)
    noise <- l == G + 1
    if(!noise) {
        if(l != G)
            stop("pro improperly specified")
        K <- G
        Vinv <- NULL
    } else {
        K <- G + 1
        Vinv <- parameters$Vinv
        if(is.null(Vinv) || Vinv <= 0)
            Vinv <- hypvol(data, reciprocal = TRUE)
    }
    if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
           any(is.null(parameters[c("pro", "mean", "variance")]))) {
        WARNING <- "parameters are missing"
        if(warn) warning(WARNING)
        z <- matrix(NA,n,K)
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        return(structure(list(modelName = "EVE", n=n, d=p, G=G, z=z,
                              parameters=parameters, loglik=NA),
                         WARNING = WARNING, returnCode = 9))
    }
    if (is.null(parameters$variance$scale) ||
            is.null(parameters$variance$shape) ||
            is.null(parameters$variance$orientation))
        stop("variance parameters are missing")
    #
    # MICHAEL from here-------------------------------------------------------
    #
    temp <- .Fortran( "eseve",
                      x = as.double(data),
                      z = double(n*K),
                      n = as.integer(n),
                      p = as.integer(p),
                      G = as.integer(G),
                      Gnoise = as.integer(K),
                      mu = as.double(mu),
                      O =  as.double( t(O) ),
                      scale = as.double(scale),
                      shape = as.double(shape),
                      pro = as.double(pro),
                      Vinv = as.double( if (is.null(Vinv)) -1 else Vinv ),
                      loglik = double(1),
                      eps = as.double(.Machine$double.eps),
                      package = "mclust")
    #
    loglik <- temp$loglik
    z <- matrix(temp$z, n,K)
    WARNING <- NULL
    if(loglik > signif(.Machine$double.xmax, 6)) {
        WARNING <- "cannot compute E-step"
        if(warn) warning(WARNING)
        z[] <- loglik <- NA
        ret <- -1
    }
    else ret <- 0
    dimnames(z) <- list(dimnames(data)[[1]],NULL)
    structure(list(modelName = "EVE", n = n, d = p, G = G,
                   z = z, parameters = parameters, loglik = loglik),
              WARNING = WARNING, returnCode = ret)
    
}

####
cdensEVE <- function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
    if(is.null(warn)) warn <- mclust.options("warn")
    dimdat <- dim(data)
    if(is.null(dimdat) || length(dimdat) != 2)
        stop("data must be a matrix")
    data <- as.matrix(data)
    n <- nrow(data)
    p <- ncol(data)
    pro <- parameters$pro
    pro <- pro/sum(pro)
    l <- length(pro)
    mu <- as.matrix(parameters$mean)
    scale <- parameters$variance$scale
    shape <- parameters$variance$shape
    O <- parameters$variance$orientation
    p <- ncol(data)
    G <- ncol(mu)
    if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
           any(is.null(parameters[c("pro", "mean", "variance")]))) {
        WARNING <- "parameters are missing"
        if(warn) warning(WARNING)
        z <- matrix(NA,n,G)
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        return(structure(z, logarithm = logarithm, modelName = "EVE",
                         WARNING = WARNING, returnCode = 9))
    }
    if (is.null(parameters$variance$scale) ||
            is.null(parameters$variance$shape) ||
            is.null(parameters$variance$orientation))
        stop("variance parameters are missing")
    #
    # MICHAEL from here-------------------------------------------------------
    #
    temp <- .Fortran( "eseve",
                      x = as.double(data),
                      z = double(n*G),
                      n = as.integer(n),
                      p = as.integer(p),
                      G = as.integer(G),
                      Gnoise = as.integer(G),
                      mu = as.double(mu),
                      O =  as.double( t(O) ),
                      scale = as.double(scale),
                      shape = as.double(shape),
                      pro = as.double(-1),
                      Vinv = as.double(-1),
                      loglik = double(1),
                      eps = as.double(.Machine$double.eps),
                      package = "mclust")
    #
    loglik <- temp$loglik
    z <- matrix(temp$z, n,G)
    WARNING <- NULL
    if(loglik > signif(.Machine$double.xmax, 6)) {
        WARNING <- "cannot compute E-step"
        if(warn) warning(WARNING)
        z[] <- NA
        ret <- -1
    } else {
        if (!logarithm) z <- exp(z)
        ret <- 0
    }
    dimnames(z) <- list(dimnames(data)[[1]],NULL)
    structure(z, logarithm = logarithm, modelName = "EVE",
              WARNING = WARNING, returnCode = ret)
}

###
simEVE <- function(parameters, n, seed = NULL, ...) 
{
    if (!is.null(seed)) 
        set.seed(seed)
    mu <- as.matrix(parameters$mean)
    d <- nrow(mu)
    G <- ncol(mu)
    if (any(is.na(parameters[c("mean", "variance")])) || 
            any(is.null(parameters[c("mean", "variance")]))) {
        warn <- "parameters are missing"
        warning("parameters are missing")
        return(structure(matrix(as.double(NA), n, d + 1), modelName = "EVE"))
    }
    pro <- parameters$pro
    if (is.null(pro)) 
        pro <- rep(1/G, G)
    clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
    ctabel <- tabulate(clabels, nbins = G)
    x <- matrix(0, n, d)
    rtshape <- sqrt(parameters$variance$shape)
    if (dim(rtshape)[1] != d | dim(rtshape)[2] != G) 
        stop("shape incompatible with mean")
    rtscale <- sqrt(parameters$variance$scale)   
    for (k in 1:G) 
    {
        m <- ctabel[k]
        sss <- rtscale * rtshape[,k]
        cholSigma <- t(parameters$variance$orientation) * sss
        x[clabels == k, ] <- sweep( matrix(rnorm(m*d), nrow = m, ncol = d) %*% cholSigma, 
                                    MARGIN = 2, STATS = mu[,k], FUN = "+" )
    }
    dimnames(x) <- list(NULL, 1:d)
    structure(cbind(group = clabels, x), modelName = "EVE")
}


##############################################################################
###                               VVE model                               ####
##############################################################################
emVVE <- function(data, parameters, prior = NULL, control = emControl(),
                  warn = NULL, ...)
{
    z <- estepVVE(data, parameters = parameters, warn = warn)$z
    meVVE(data, z = z, prior = prior, control = control,
          Vinv = parameters$Vinv, warn = warn)
}

####
meVVE <- function(data, z, prior = NULL, control = emControl(),
                  Vinv = NULL, warn = NULL, ...)
{
    if(is.null(warn)) warn <- mclust.options("warn")
    dimdat <- dim(data)
    oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
    if(oneD || length(dimdat) != 2)
        stop("data should in the form of a matrix")
    data <- as.matrix(data)
    n <- nrow(data)
    p <- ncol(data)
    z <- as.matrix(z)
    dimz <- dim(z)
    if(dimz[1] != n)
        stop("data and z should have the same row dimension")
    K <- dimz[2]
    if (!is.null(Vinv)) {
        G <- K - 1
        if(Vinv <= 0) Vinv <- hypvol(data, reciprocal = TRUE)
    } else G <- K
    if(all(is.na(z))) {
        WARNING <- "z is missing"
        if(warn) warning(WARNING)
        variance <- list(modelName = "VVE", d = p, G = G,
                         scale=rep(NA,G), shape=rep(NA,p), orientation=array(NA,c(p,p,G)))
        parameters <- list(pro=rep(NA,G), mean=matrix(NA,p,G),
                           variance=variance)
        return(structure(list(modelName="VVE", prior=prior, n=n, d=p,
                              G=G, z=z, parameters=parameters,
                              control=control, loglik=NA),
                         WARNING = WARNING, returnCode = 9))
    }
    if(any(is.na(z)) || any(z < 0) || any(z > 1))
        stop("improper specification of z")
    lwork <- max(3 * min(n, p) + max(n, p), 5 * min(n, p), p + G)
    storage.mode(z) <- "double"

    # MICHAEL from here-------------------------------------------------------
    #
    # without prior specification
    if(is.null(prior)) {
        temp <- .Fortran("mevve",
                         x = as.double(data),
                         z = as.double(z),
                         n = as.integer(n),
                         p = as.integer(p),
                         G = as.integer(G),
                         Gnoise = as.integer(K),
                         mu = double(p*G),
                         O = as.double( diag(p) ),
                         U = double(p*p*G),
                         scale = as.double( rep(1, G) ),
                         shape = as.double( matrix(1, p,G) ),
                         pro = double(K),
                         Vinv = as.double( if (is.null(Vinv)) -1 else Vinv ),
                         loglik = double(1),
                         eqpro = as.logical(control$equalPro),
                         itmaxin = as.integer(control$itmax[2]),
                         tolin = as.double(control$tol[2]),
                         itmaxout = as.integer(control$itmax[1]),
                         tolout = as.double(control$tol[1]),
                         eps = as.double(control$eps),
                         niterin = integer(1),
                         errin = double(1),
                         niterout = integer(1),
                         errout = double(1),
                         lwork = as.integer(lwork),
                         info = as.integer(0),
                         package = "mclust")
        #
    } 
    else {
        # with prior
        priorParams <- do.call(prior$functionName,
                               c(list(data = data, G = G, modelName = "VVE"),
                                 prior[names(prior) != "functionName"]))
        # temp <- .Fortran("mevvep", ...)
        temp <- list(x = data, z = z, n = n, p = p, G = G,
                     mu = double(p*G), O = double(p*p), U = double(p*p*G),
                     scale = as.double(rep(1, G)), shape = double(p*G), 
                     pro = double(G), loglik = NA, 
                     eqpro = as.logical(control$equalPro),
                     itmaxin = as.integer(control$itmax[2]),
                     tolin = as.double(control$tol[2]),
                     itmaxout = as.integer(control$itmax[1]),
                     tolout = as.double(control$tol[1]),
                     eps = as.double(control$eps),
                     niterin = integer(1), errin = double(1),
                     niterout = integer(1), errout = double(1),
                     lwork = as.integer(lwork), info = FALSE)
        WARNING <- "VVE model is not available with prior"
        if(warn) warning(WARNING)
        temp <- structure(temp, info = NA, WARNING = WARNING, returnCode = -1)
        return(temp)
    }
    z <- matrix(temp$z, n,K)
    niterin <- temp$niterin
    errin <- temp$errin
    niterout <- temp$niterout
    errout <- temp$errout
    loglik <- temp$loglik
    lapackSVDinfo <- temp$info
    mu <- matrix(temp$mu, p,G)
    dimnames(mu) <- list(NULL, as.character(1:G))
    scale <- temp$scale
    shape <- matrix(temp$shape, p,G)
    O <- t( matrix(temp$O, p,p) )
    pro <- temp$pro
    if( !is.finite(loglik) | 
        any(scale > signif(.Machine$double.xmax, 6)) |
        any(shape > signif(.Machine$double.xmax, 6)) |
        any(O > signif(.Machine$double.xmax, 6)) |
        any(is.nan(scale)) | any(is.nan(shape)) | any(is.nan(O)) )
      { loglik <- .Machine$double.xmax }
    WARNING <- NULL
    #
    if(lapackSVDinfo) {
        if(lapackSVDinfo > 0) {
            WARNING <- "LAPACK DSYEV or DGESVD fails to converge"
        }
        else {
            WARNING <- "input error for LAPACK DSYEV or DGESVD"
        }
        if(warn) warning(WARNING)
        O[] <- shape[] <- scale[] <- NA
        mu[] <- pro[] <- z[] <- loglik <- NA
        sigma <- array(NA, c(p, p, G))
        ret <- -9
        #
    } else if(loglik > signif(.Machine$double.xmax, 6)) {
        WARNING <- "singular covariance"
        if(warn) warning(WARNING)
        O[] <- shape[] <- scale[] <- NA
        mu[] <- pro[] <- z[] <- loglik <- NA
        sigma <- array(NA, c(p, p, G))
        ret <- -1
    } else if(loglik <  - signif(.Machine$double.xmax, 6)) {
        if(control$equalPro) {
            WARNING <- "z column sum fell below threshold"
            if(warn) warning(WARNING)
        } else {
            WARNING <- "mixing proportion fell below threshold"
            if(warn) warning(WARNING)
        }
        mu[] <- pro[] <- z[] <- loglik <- NA
        sigma <- array(NA, c(p, p, G))
        ret <- if(control$equalPro) -2 else -3
    } else {
        sigma <- array(NA, c(p,p,G))
        for ( g in 1:G ) sigma[,,g] <- scale[g] * O %*% diag(shape[,g]) %*% t(O)
        if(niterin >= control$itmax[2]) {
            WARNING <- "inner iteration limit reached"
            if(warn) warning(WARNING)
            niterin <-  - niterin
            ret <- 2
        } else if(niterout >= control$itmax[1]) {
            WARNING <- "iteration limit reached"
            if(warn) warning(WARNING)
            niterout <-  - niterout
            ret <- 1
        } else ret <- 0
    }
    info <- structure(c(niterout = niterout, errout = errout),
                      inner = c(niterin = niterin, errin = errin))
    # info <- structure(c(iterations = its, error = err),
    #                   inner = c(iterations = inner, error = inerr))
    dimnames(z) <- list(dimnames(data)[[1]],NULL)
    dimnames(mu) <- list(dimnames(data)[[2]], NULL)
    dimnames(sigma) <- dimnames(O) <- list(dimnames(data)[[2]], dimnames(data)[[2]])
    ##  Sigma = scale * O %*% diag(shape) %*% t(O)
    variance <- list(modelName = "VVE", d = p, G = G, sigma = sigma,
                     scale = scale, shape = shape, orientation = O)
    parameters <- list(pro=pro, mean=mu, variance=variance, Vinv=Vinv)
    structure(list(modelName = "VVE", prior = prior, n = n, d = p, G = G,
                   z = z, parameters = parameters, control = control,
                   loglik = loglik),
              info = info, WARNING = WARNING, returnCode = ret)
}


####
mstepVVE <- function(data, z, prior = NULL, warn = NULL, control = NULL, ...)
{
    if (is.null(warn)) warn <- mclust.options("warn")
    dimdat <- dim(data)
    oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
    if(oneD || length(dimdat) != 2)
        stop("data should be a matrix or a vector")
    data <- as.matrix(data)
    n <- nrow(data)
    p <- ncol(data)
    ##
    z <- as.matrix(z)
    dimz <- dim(z)
    if(dimz[1] != n)
        stop("row dimension of z should equal data length")
    G <- dimz[2]
    if(all(is.na(z))) {
        WARNING <- "z is missing"
        if(warn) warning(WARNING)
        variance <- list(modelName = "VVE", d = p, G = G,
                         scale = rep(NA,G), shape = rep(NA,p), orientation = array(NA,c(p,p,G)))
        parameters <- list(pro=rep(NA,G), mean=matrix(NA,p,G),
                           variance=variance)
        return(structure(list(modelName="VVE", prior=prior, n=n, d=p,
                              G=G, z=z, parameters=parameters,
                              control=control, loglik=NA),
                         WARNING = WARNING, returnCode = 9))
        
        WARNING <- "z is missing"
        if(warn) warning(WARNING)
        return(structure(list(
            n = n, d = p, G = G, mu = matrix(NA,p, G), sigma = array(NA, c(p, p, G)),
            decomp = list(d = p, G = G, scale = rep(NA, G), shape = rep(NA, p),
                          orientation = array(NA, c(p, p, G))),
            pro = rep(NA,G), modelName = "VVE", prior = prior), WARNING = WARNING))
    }
    if(any(is.na(z)) || any(z < 0) || any(z > 1)) stop("improper specification of z")
    if (is.null(control)) control <- emControl()
    itmax <- if(length(control$itmax) == 1) control$itmax else control$itmax[2]
    tol <- if(length(control$tol) == 1) control$tol else control$tol[2]
    lwork <- max(3 * min(n, p) + max(n, p), 5 * min(n, p), p + G)
    
    #
    # MICHAEL from here-------------------------------------------------------
    #
    # without prior specification
    if(is.null(prior)) {
        temp <- .Fortran("msvve",
                         x = as.double(data), 
                         z = as.double(z),
                         n = as.integer(n),
                         p = as.integer(p),
                         G = as.integer(G),
                         mu = double(p*G),
                         U = double(p*p*G),
                         O = as.double( diag(p) ),
                         scale = as.double( rep(1, G) ),
                         shape = as.double( matrix(1, p,G) ),
                         pro = double(G),
                         lwork = as.integer(lwork),
                         info = as.integer(0),
                         itmax = as.integer(itmax),
                         tol = as.double(tol),
                         niterin = integer(1),
                         errin = double(1),
                         eps = as.double(.Machine$double.eps),
                         package = "mclust")
    } else {
        # with prior
        priorParams <- do.call(prior$functionName,
                               c(list(data = data, G = G, modelName = "VVE"),
                                 prior[names(prior) != "functionName"]))
        #
        # temp <- .Fortran("msvvep", ...)
        temp <- list(x = data, z = z, n = n, p = p, G = G,
                     mu = double(p*G),  U = double(p*p*G), O = double(p*p),
                     scale = double(1), pro = double(G), shape = double(p*G),
                     lwork = as.integer(lwork), info = FALSE,
                     itmax = as.integer(itmax), tol = as.double(tol),
                     niterin = integer(1), errin = double(1),
                     eps = as.double(.Machine$double.eps))
        WARNING <- "VVE model is not available with prior"
        if(warn) warning(WARNING)
    }
    lapackSVDinfo <- temp$info
    errin <- temp$errin
    niterin <- temp$niterin
    mu <- matrix(temp$mu, p,G)
    dimnames(mu) <- list(NULL, as.character(1:G))
    O <- t( matrix(temp$O, p,p) )
    shape <- matrix(temp$shape, p,G)
    scale <- temp$scale
    pro <- temp$pro
    WARNING <- NULL
    #
    if(lapackSVDinfo) {
        if(lapackSVDinfo > 0) {
            WARNING <- "LAPACK DSYEV or DGESVD fails to converge"
        } else {
            WARNING <- "input error for LAPACK DSYEV or DGESVD"
        }
        if(warn) warning(WARNING)
        O[] <- shape[] <- scale[] <- NA
        sigma <- array(NA, c(p, p, G))
        ret <- -9
        #
    } else if(any(c(scale, shape) > signif(.Machine$double.xmax, 6))) {
        WARNING <- "cannot compute M-step"
        if(warn) warning(WARNING)
        mu[] <- pro[] <- O[] <- shape[] <- scale[] <- NA
        sigma <- array(NA, c(p, p, G))
        ret <- -1
    } else {
        # sigma <- array( apply(shape, 2, function(sh) O%*%diag(sh)%*%t(O)), c(p,p,G) )
        sigma <- array(NA, c(p,p,G))
        for ( g in 1:G ) sigma[,,g] <- scale[g] * O %*% diag(shape[,g]) %*% t(O)
        if(niterin >= itmax) {
            WARNING <- "inner iteration limit reached"
            if(warn) warning(WARNING)
            niterin <-  - niterin
        }
        ret <- 2
    }
    info <- c(iteration = niterin, error = errin)
    dimnames(z) <- list(dimnames(data)[[1]], NULL)
    dimnames(mu) <- list(dimnames(data)[[2]], NULL)
    dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]])
    dimnames(O) <- list(dimnames(data)[[2]], dimnames(data)[[2]])
    variance <- list(modelName = "VVE", d = p, G = G, sigma = sigma,
                     scale = scale, shape = shape, orientation = O)
    parameters <- list(pro=pro, mean=mu, variance=variance)
    structure(list(modelName = "VVE", prior = prior, n = n, d = p, G = G,
                   z = z, parameters = parameters, control = control),
              info = info, WARNING = WARNING, returnCode = ret)
}


###
estepVVE <- function(data, parameters, warn = NULL, ...)
{
    if (is.null(warn)) warn <- mclust.options("warn")
    dimdat <- dim(data)
    if(is.null(dimdat) || length(dimdat) != 2)
        stop("data must be a matrix")
    data <- as.matrix(data)
    n <- nrow(data)
    p <- ncol(data)
    pro <- parameters$pro
    pro <- pro/sum(pro)
    l <- length(pro)
    mu <- as.matrix(parameters$mean)
    scale <- parameters$variance$scale
    shape <- parameters$variance$shape
    O <- parameters$variance$orientation
    G <- ncol(mu)
    noise <- l == G + 1
    if(!noise) {
        if(l != G)
            stop("pro improperly specified")
        K <- G
        Vinv <- NULL
    } else {
        K <- G + 1
        Vinv <- parameters$Vinv
        if(is.null(Vinv) || Vinv <= 0)
            Vinv <- hypvol(data, reciprocal = TRUE)
    }
    if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
           any(is.null(parameters[c("pro", "mean", "variance")]))) {
        WARNING <- "parameters are missing"
        if(warn) warning(WARNING)
        z <- matrix(NA,n,K)
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        return(structure(list(modelName = "VVE", n=n, d=p, G=G, z=z,
                              parameters=parameters, loglik=NA),
                         WARNING = WARNING, returnCode = 9))
    }
    if (is.null(parameters$variance$scale) ||
            is.null(parameters$variance$shape) ||
            is.null(parameters$variance$orientation))
        stop("variance parameters are missing")
    #
    # MICHAEL from here-------------------------------------------------------
    #
    temp <- .Fortran( "esvve",
                      x = as.double(data),
                      z = double(n*K),
                      n = as.integer(n),
                      p = as.integer(p),
                      G = as.integer(G),
                      Gnoise = as.integer(K),
                      mu = as.double(mu),
                      O =  as.double( t(O) ),
                      scale = as.double(scale),
                      shape = as.double(shape),
                      pro = as.double(pro),
                      Vinv = as.double( if (is.null(Vinv)) -1 else Vinv ),
                      loglik = double(1),
                      eps = as.double(.Machine$double.eps),
                      package = "mclust")
    #
    loglik <- temp$loglik
    z <- matrix(temp$z, n,K)
    WARNING <- NULL
    if(loglik > signif(.Machine$double.xmax, 6)) {
        WARNING <- "cannot compute E-step"
        if(warn) warning(WARNING)
        z[] <- loglik <- NA
        ret <- -1
    }
    else ret <- 0
    dimnames(z) <- list(dimnames(data)[[1]],NULL)
    structure(list(modelName = "VVE", n = n, d = p, G = G,
                   z = z, parameters = parameters, loglik = loglik),
              WARNING = WARNING, returnCode = ret)
    
}

####
cdensVVE <- function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
    if(is.null(warn)) warn <- mclust.options("warn")
    dimdat <- dim(data)
    if(is.null(dimdat) || length(dimdat) != 2)
        stop("data must be a matrix")
    data <- as.matrix(data)
    n <- nrow(data)
    p <- ncol(data)
    pro <- parameters$pro
    pro <- pro/sum(pro)
    l <- length(pro)
    mu <- as.matrix(parameters$mean)
    scale <- parameters$variance$scale
    shape <- parameters$variance$shape
    O <- parameters$variance$orientation
    p <- ncol(data)
    G <- ncol(mu)
    if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
           any(is.null(parameters[c("pro", "mean", "variance")]))) {
        WARNING <- "parameters are missing"
        if(warn) warning(WARNING)
        z <- matrix(NA,n,G)
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        return(structure(z, logarithm = logarithm, modelName = "VVE",
                         WARNING = WARNING, returnCode = 9))
    }
    if (is.null(parameters$variance$scale) ||
            is.null(parameters$variance$shape) ||
            is.null(parameters$variance$orientation))
        stop("variance parameters are missing")
    #
    # MICHAEL from here-------------------------------------------------------
    #
    temp <- .Fortran( "esvve",
                      x = as.double(data),
                      z = double(n*G),
                      n = as.integer(n),
                      p = as.integer(p),
                      G = as.integer(G),
                      Gnoise = as.integer(G),
                      mu = as.double(mu),
                      O =  as.double( t(O) ),
                      scale = as.double(scale),
                      shape = as.double(shape),
                      pro = as.double(-1),
                      Vinv = as.double(-1),
                      loglik = double(1),
                      eps = as.double(.Machine$double.eps),
                      package = "mclust")
    #
    loglik <- temp$loglik
    z <- matrix(temp$z, n,G)
    WARNING <- NULL
    if(loglik > signif(.Machine$double.xmax, 6)) {
        WARNING <- "cannot compute E-step"
        if(warn) warning(WARNING)
        z[] <- NA
        ret <- -1
    } else {
        if (!logarithm) z <- exp(z)
        ret <- 0
    }
    dimnames(z) <- list(dimnames(data)[[1]],NULL)
    structure(z, logarithm = logarithm, modelName = "VVE",
              WARNING = WARNING, returnCode = ret)
}

###
simVVE <- function(parameters, n, seed = NULL, ...) 
{
    if (!is.null(seed)) 
        set.seed(seed)
    mu <- as.matrix(parameters$mean)
    d <- nrow(mu)
    G <- ncol(mu)
    if (any(is.na(parameters[c("mean", "variance")])) || 
            any(is.null(parameters[c("mean", "variance")]))) {
        warn <- "parameters are missing"
        warning("parameters are missing")
        return(structure(matrix(as.double(NA), n, d + 1), modelName = "VVE"))
    }
    pro <- parameters$pro
    if (is.null(pro)) 
        pro <- rep(1/G, G)
    clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
    ctabel <- tabulate(clabels, nbins = G)
    x <- matrix(0, n, d)
    rtshape <- sqrt(parameters$variance$shape)
    if (dim(rtshape)[1] != d | dim(rtshape)[2] != G) 
        stop("shape incompatible with mean")
    rtscale <- sqrt(parameters$variance$scale)
    if (length(rtscale) != G) 
        stop("scale incompatible with mean")
    for (k in 1:G) 
    {
        m <- ctabel[k]
        sss <- rtscale[k] * rtshape[,k]
        cholSigma <- t(parameters$variance$orientation) * sss
        x[clabels == k, ] <- sweep( matrix(rnorm(m*d), nrow = m, ncol = d) %*% cholSigma, 
                                    MARGIN = 2, STATS = mu[,k], FUN = "+" )
    }
    dimnames(x) <- list(NULL, 1:d)
    structure(cbind(group = clabels, x), modelName = "VVE")
}
