####################################################
### Authors: Simone Padoan and Moreno Bevilacqua.
### Emails: simone.padoan@unibocconi.it,
### moreno.bevilacqua@uv.cl
### Institutions: Department of Decision Sciences,
### University Bocconi of Milan and
### Departamento de Estadistica
### Universidad de Valparaiso
### File name: Distribution.r
### Description:
### This file contains a set of procedures
### that define several distributions and
### related functions.
### Last change: 28/03/2013.
####################################################


### Procedures are in alphabetical order.

FitGev <- function(data, method='Nelder-Mead', start, varest=FALSE)
  {
    param <- NULL
    varcov <- NULL
    stderr <- NULL
    numdata <- length(data)

    if(missing(data) || numdata == 0 || !is.numeric(data))
      stop('insert a numeric vector of data')

    ### Starting values:
    momest <- MomEst(data, numdata)
    location <- momest$location
    scale <- momest$scale
    shape <- momest$shape
    if(shape >= 0) shape <- .1
    else shape <- -.1

    param <- list(location=location, scale=scale, shape=shape)
    parnames <- c('location', 'scale', 'shape')

    if(!missing(start))
      {
        if(!is.list(start))
          stop('start need to be named list of starting values')

        for(i in 1 : length(start))
          {
            namestart <- names(start[i])
            if(!any(namestart == parnames))
              stop('insert location, scale and shape as a named list')

            param[namestart] <- start[namestart]

          }

        if(param$scale < 0)
          param$scale <- scale
      }

    param <- unlist(param)

    ### fit the log-likelihood of the GEV model:
    fitted <- optim(param, GevLogLik, data=data, method=method, numdata=numdata,
                    control=list(fnscale=-1, reltol=1e-14, maxit=1e8), hessian=varest)

    if(varest)
      {
        varcov <- - try(solve(fitted$hessian), silent = TRUE)
        if(!is.matrix(varcov))
          {
            varcov <- 'none'
            stderr <- 'none'
          }
        else
          {
            stderr <- diag(varcov)
            if(any(stderr < 0))
              stderr <- 'none'
            else
              stderr <- sqrt(stderr)
          }
      }

    return(list(param=fitted$par, varcov=varcov, stderr=stderr))
  }

Dist2Dist <- function(data, from='Gev', to='sFrechet', loc=NULL, scale=NULL, shape=NULL)
  {
    Dist2Dist <- NULL

    if(missing(data) || !is.numeric(data))
      stop('insert a numeric vector of data')

    dimdata <- dim(data)
    if(is.null(dimdata))
      {
        numdata <- length(data)
        if(numdata==0)
          stop('insert a numeric vector of data')
        numcoord <- 1
        dimdata <- c(numdata, numcoord)
        dim(data) <- dimdata
        numparam <- numdata
      }
    else
      {
        if(length(dimdata)==2)
          {
            numdata <- dimdata[1]
            numcoord <- dimdata[2]
          }
        if(length(dimdata)==3)
          {
            numdata <- dimdata[3]
            numcoord <- prod(dimdata[1:2])
          }
        numparam <- numcoord
      }
    if(is.null(loc)) loc <- double(numcoord)
    else if(!is.numeric(loc) || !length(loc)==numparam)
      stop('insert a numeric vector of locations\n')
    if(is.null(scale)) scale <- double(numcoord) + 1
    else if(!is.numeric(scale) || !length(scale)==numparam)
      stop('insert a numeric vector of scales\n')
    if(is.null(shape)) shape <- double(numcoord)
    else if(!is.numeric(shape) || !length(shape)==numparam)
      stop('insert a numeric vector of shapes\n')
    # Define the fitting procedure to compute the ML estimates:
    Mle <- function(x)
      return(as.double(FitGev(x)$param))
    # Transform from GEV to other distributions:
    if(from=='Gev')
      {
        type <- switch(to,
                       Uniform=0,
                       sFrechet=1,
                       sGumbel=2,
                       sWeibull=3,
                       Gev=4)
        if(is.null(type))
          stop('insert a valid distribution name\n')
        param <- apply(data, 2, Mle) # ML estimates
        # Margins transformation:
        Dist2Dist <- .C('Dist2Dist', as.double(data), as.double(param[1,]), as.double(param[2,]),
                        as.double(param[3,]), as.integer(numdata), as.integer(numcoord),
                        as.double(loc), as.double(scale), as.double(shape), as.integer(type),
                        res=double(numdata * numcoord), PACKAGE='CompRandFld', DUP=TRUE, NAOK=TRUE)$res
      }
    # Transform from standard Frechet to the GEV distribution:
    if(from=='sFrechet')
      {
        type <- switch(to, Gev=5)
        if(is.null(type))
          stop('insert a valid distribution name\n')
        param <- rep(1, 3)
        # Margins transformation:
        Dist2Dist <- .C('Dist2Dist', as.double(data), as.double(param[1]), as.double(param[2]),
                        as.double(param[3]), as.integer(numdata), as.integer(numcoord),
                        as.double(loc), as.double(scale), as.double(shape), as.integer(type),
                        res=double(numdata * numcoord), PACKAGE='CompRandFld', DUP=TRUE, NAOK=TRUE)$res
      }
    # Transform from standard Gumbel to the GEV distribution:
    if(from=='sGumbel')
      {
        type <- switch(to, Gev=6)
        if(is.null(type))
          stop('insert a valid distribution name\n')
        param <- rep(1, 3)
        # Margins transformation:
        Dist2Dist <- .C('Dist2Dist', as.double(data), as.double(param[1]), as.double(param[2]),
                        as.double(param[3]), as.integer(numdata), as.integer(numcoord),
                        as.double(loc), as.double(scale), as.double(shape), as.integer(type),
                        res=double(numdata * numcoord), PACKAGE='CompRandFld', DUP=TRUE, NAOK=TRUE)$res
      }
    if(from=='sWeibull')
      {
        type <- switch(to, Gev=7)
        if(is.null(type))
          stop('insert a valid distribution name\n')
        param <- rep(1, 3)
        # Margins transformation:
        Dist2Dist <- .C('Dist2Dist', as.double(data), as.double(param[1]), as.double(param[2]),
                        as.double(param[3]), as.integer(numdata), as.integer(numcoord),
                        as.double(loc), as.double(scale), as.double(shape), as.integer(type),
                        res=double(numdata * numcoord), PACKAGE='CompRandFld', DUP=TRUE, NAOK=TRUE)$res
      }

    dim(Dist2Dist) <- dimdata

    return(Dist2Dist)
  }

GevLogLik <- function(data, numdata, param)
  {
    ### Compute the log-likelihood:
    result <- .C('GevLogLik', as.double(data), as.integer(numdata), as.double(param),
                 res=double(1), PACKAGE='CompRandFld', DUP=TRUE, NAOK=TRUE)$res
    return(result)
  }

MomEst <- function(data, n)
  {
    # Scale estimate:
    scale <- sqrt(6 * var(data)) / pi
    # Location estimate:
    location <- mean(data) - 0.58 * scale

    k <- round(.5 * n)
    mind <- min(data)
    if(mind < 0) data <- data - mind

    data <- sort(data)
    delta <- log(data[n-0:(k-1)] / data[n-k])
    M1 <- sum(delta) / k
    M2 <- sum(delta^2) / k
    # Shape estimate:
    shape <- M1 + 1 - .5 * (M2 / (M2 - M1^2))

    return(list(location=location, scale=scale, shape=shape))
  }
