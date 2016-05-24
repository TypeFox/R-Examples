## --------------------------------------------------------------------------
##
## Table of x-errors
##
## --------------------------------------------------------------------------

xerror <- function (n, aqdist, qdist, ..., trunc=NULL, udomain=c(0,1),
                    res=1000, kind=c("abs","rel"), tails=FALSE, plot=FALSE)

  ## ------------------------------------------------------------------------
  ## Create table of x-errors for numerical inversion method.
  ## ------------------------------------------------------------------------
  ## n      : Size of random sample 
  ## aqdist : Approximate inverse distribution function (quantile function):
  ##             a function or an object of class "unuran"
  ## qdist  : (Exact) quatile function of distribution
  ## ....   : Parameters of distribution
  ## trunc  : boundaries of truncated domain of distribution
  ## udomain: domain for u
  ## res    : Resolution of table (number of intervals in [0,1] for which 
  ##          quantiles a sample of x-errors are computed and stored) 
  ## kind   : kind of x-error: absolute, relative
  ## tails  : if TRUE, then the tail regions are treated more accurately
  ## plot   : if TRUE, then plot a u-errors
  ## ------------------------------------------------------------------------
  ## Return:
  ##   list of size 'res' that contains 'n', 'res', domain, kind and
  ##   the quantiles of x-errors:
  ##     min, lqr (lower quartile), med, uqr (upper quartile), max 
  ## ------------------------------------------------------------------------
{
  ## sample size
  if (missing(n) || !is.numeric(n) || n<1 || n!=round(n))
    stop ("Invalid argument 'n'.")

  ## resolution
  if (!is.numeric(res) || res<1 || res!=round(res))
    stop ("Invalid argument 'res'.")

  ## kind of x-error
  kind <- match.arg(kind)

  kind.name <- switch(kind,
                      "abs" = "absolute x-error",
                      "rel" = "relative x-error",
                      stop ("invalid 'kind'") )

  ## approximate inverse distribution function (quantile function)
  if( missing(aqdist) )
    stop ("Argument 'aqdist' missing.")

  if (is.function(aqdist)) {
    ## R function
    myaqdist <- aqdist
  }
  else if (is(aqdist,"unuran")) {
    ## "unuran" object
    ## Remark: We assume that package 'Runuran' is already loaded
    ##    because it is required to create an object of class "unuran".
    ## However, we need an object that contains a
    ## univariate continuous distribution.
    if (unuran.distr.class(aqdist) != "cont")
      stop ("Argument 'aqdist' is object of class 'unuran' of invalid distribution type.")
    ## Moreover, we need an inversion method
    if (!unuran.is.inversion(aqdist))
      stop ("Argument 'aqdist' is invalid UNU.RAN object: inversion method required.")

    myaqdist <- function(u) { uq(unr=aqdist, u) }
  }
  else {
    stop ("Argument 'aqdist' invalid.")
  }
  
  ## (exact) quatile function of distribution
  if( missing(qdist) || !is.function(qdist))
    stop ("Argument 'qdist' missing or invalid.")

  ## handle truncated domain (for target distribution)
  if (! is.null(trunc) ) {
    if (! (length(trunc)==2 && trunc[1]<trunc[2]))
      stop ("Argument 'trunc' invalid.")

    CDFmin <- invqdist(trunc[1], qdist, ...)
    CDFmax <- invqdist(trunc[2], qdist, ...)

    tmpq <- qdist
    qdist <- function(x) { tmpq(x * (CDFmax - CDFmin) + CDFmin, ...) }
  }
  
  ## domain for quantile function
  umin <- max(0,udomain[1])
  umax <- min(1,udomain[2])
  if( umin>=umax) 
    stop ("Invalid argument 'udomain'.")

  ## we do not treat tails more accurately
  ## when a smaller domain is given
  if(!isTRUE(all.equal(c(umin,umax),c(0,1))))
    tails <- FALSE
  
  ## interval boundaries
  length <- (umax-umin) / res
  uiv <- umin + length * (0:res)
  
  ## samplesize for each interval
  k <- round(n / res)
  if( k<1 ) stop ("Invalid arguments 'n' < 'res'.")
  n <- k * res

  ## arrays for storing u-errors
  xe.min <- numeric(res)
  xe.lqr <- numeric(res)
  xe.med <- numeric(res)
  xe.uqr <- numeric(res)
  xe.max <- numeric(res)
  xe.mad <- numeric(res)
  xe.mse <- numeric(res)

  ## loop over intervals
  for (i in 1:res) {

    ## u values
    if (isTRUE(tails) && ( isTRUE(all.equal(i,1)) || isTRUE(all.equal(i,res)) )) {
      DBL.EPSILON <- 2^(-52)
      size <- round(n/2)
      if (isTRUE(all.equal(i,1))) {
        ivmin <- DBL.EPSILON
        ivmax <- 1/res
      }
      else {
        ivmin <- 1 - 1/res
        ivmax <- 1 - DBL.EPSILON
      }
      u <- ivmin + (0:(size-1)) * (ivmax - ivmin)/(size-1)
    }
    else {
      u <- uiv[i] + length * ((1:k) - 0.5)/k
    }

    ## compute x-error
    if (kind == "abs") {
      Fu <- qdist(u,...)
      xerr <- abs(Fu - myaqdist(u))
    }
    else if (kind == "rel") {
      Fu <- qdist(u,...)
      xerr <- abs(Fu - myaqdist(u)) / abs(Fu)
    }
    else {
      stop ("invalid 'kind'")
    }
    
    ## quantiles for error
    xdata <- quantile(xerr,c(0,0.25,0.5,0.75,1))

    xe.min[i] <- xdata[[1]]
    xe.lqr[i] <- xdata[[2]]
    xe.med[i] <- xdata[[3]]
    xe.uqr[i] <- xdata[[4]]
    xe.max[i] <- xdata[[5]]

    ## MAD and MSE
    xe.mad[i] <- sum(xerr) / length(xerr)
    xe.mse[i] <- sum(xerr^2) / length(xerr)
  }
  
  ## return result as object of class "rvgt.ierror"
  xerror <- list(n=n, res=res, udomain=c(umin,umax), kind=kind.name,
                 min=xe.min, lqr=xe.lqr, med=xe.med, uqr=xe.uqr, max=xe.max,
                 mad=xe.mad, mse=xe.mse)
  class(xerror) <- "rvgt.ierror"
  
  ## plot x-errors
  if( isTRUE(plot) ) {
    plot(xerror)
  }

  ## return table of x-errors
  return(invisible(xerror))
}

## --------------------------------------------------------------------------
