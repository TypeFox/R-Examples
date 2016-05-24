## --------------------------------------------------------------------------
##
## RVG frequency table
##
## --------------------------------------------------------------------------

rvgt.ftable <- function (n, rep=1, rdist, qdist, pdist, ...,
                         breaks=101, trunc=NULL, exactu=FALSE, plot=FALSE)

  ## ------------------------------------------------------------------------
  ## Create RVG frequency table for random variates generator.
  ## Each row contains the frequencies for one sample of size n.
  ## ------------------------------------------------------------------------
  ## n      : Size of random sample at each repetition
  ## rep    : Number of repetitions
  ## rdist  : Random variate generator: a function or an object of class "unuran"
  ## qdist  : Quantile function of distribution
  ## pdist  : Cumulative distribution function of distribution
  ## ....   : Parameters of distribution
  ## breaks : A single number giving the number of cells of histogram; or
  ##          a vector giving the breakpoints between histogram cells
  ##          (in u-scale)
  ## trunc  : boundaries of truncated domain 
  ## exactu : Whether exact location of break points in u-scale must be used.
  ##          If FALSE, then break points are slightly moved in order of
  ##          faster runtimes (this does not effect correctness of the
  ##          frequency table.)
  ##          If TRUE this might be quite slow unless 'qdist' is given
  ##          (only if the number break points are given,
  ##          not a vector of breaks points.) 
  ## plot   : Whether to plot a histogram
  ## ------------------------------------------------------------------------
  ## return:
  ## (rep x (#breaks-1))-matrix of frequencies where each row contains
  ##   frequencies of sample of size n
  ## ------------------------------------------------------------------------
{
  ## --- constants ----------------------------------------------------------

  min.bin.width <- 1e-12   ## minimal width for bins

  ## --- variables ----------------------------------------------------------

  dtype <- NULL

  ## --- check arguments ----------------------------------------------------

  ## sample size
  if (missing(n) || !is.numeric(n) || n<100 || n!=round(n))
    stop ("Argument 'n' missing or invalid.")

  ## number of repetitions
  if (!is.numeric(rep) || rep<1 || rep!=round(rep))
    stop ("Invalid argument 'rep'.")

  ## random variate generator
  if (missing(rdist))
    stop ("Argument 'rdist' missing.")

  if (is.function(rdist)) {
    ## R function
    myrdist <- function(size) { rdist(size,...) }
  }
  else if (is(rdist,"unuran")) {
    ## "unuran" object
    ## Remark: We assume that package 'Runuran' is already loaded
    ##    because it is required to create an object of class "unuran".
    ## However, we need an object that contains a
    ## univariate continuous distribution.

    ## store type of distribution
    dtype <- unuran.distr.class(rdist)

    ## check type of distribution
    if (! (dtype == "cont" || dtype == "discr"))
      stop ("Argument 'rdist' is object of class 'unuran' of invalid distribution type.")

    ## define sampling routine
    myrdist <- function(size) { ur(unr=rdist, size) }
  }
  else {
    stop ("Argument 'rdist' invalid.")
  }
  
  ## quantile and distribution function
  if (missing(qdist)) qdist <- NULL
  if (missing(pdist)) pdist <- NULL
  if (is.null(qdist) && is.null(pdist))
    stop ("Argument 'qdist' or 'pdist' required.")

  if (!is.null(qdist) && !is.function(qdist))
    stop ("Argument 'qdist' invalid.")
  if( !is.null(pdist) && !is.function(pdist))
    stop ("Argument 'pdist' invalid.")

  ## break points
  if (!is.numeric(breaks) || length(breaks) < 1)
    stop ("Invalid argument 'breaks'.")

  ## use exact location of break points
  if (!is.logical(exactu))
    stop ("Argument 'exactu' must be boolean.")

  ## --- handle truncated domain --------------------------------------------

  if (! is.null(trunc) ) {
    if (! (length(trunc)==2 && trunc[1]<trunc[2]))
      stop ("Argument 'trunc' invalid.")

    if (! is.null(pdist)) {
      CDFmin <- pdist(trunc[1],...)
      CDFmax <- pdist(trunc[2],...)
    }
    else {
      CDFmin <- invqdist(trunc[1], qdist, ...)
      CDFmax <- invqdist(trunc[2], qdist, ...)
    }

    if (! is.null(pdist)) {
      tmpp <- pdist
      pdist <- function(x) { (tmpp(x,...) - CDFmin) / (CDFmax - CDFmin) }
    }
    if (! is.null(qdist)) {
      tmpq <- qdist
      qdist <- function(x) { tmpq(x * (CDFmax - CDFmin) + CDFmin, ...) }
    }
  }
  
  ## --- compute break points in u-scale ------------------------------------
  
  ## case: number of break points
  if (length(breaks) == 1) {
    breaks <- as.integer(breaks)
    if (breaks < 3) 
      stop (paste("Number of break points too small (less than 3):",breaks))

    ## equidistributed break points for uniform scale
    ubreaks <- (0:(breaks-1))/(breaks-1)
  }

  ## case: vector of break points (in u-scale)
  else {
    if (length(breaks) < 3) 
      stop (paste("Number of break points too small (less than 3):",length(breaks)))
    if (min(breaks)<0 || max(breaks)>1)
      stop ("break points out of range [0,1]")

    ## the break points must be sorted
    ubreaks <- sort(breaks)
    ## first and last break point must be 0 and 1, resp.
    ubreaks[1] <- 0
    ubreaks[length(ubreaks)] <- 1

    ## differences must be strictly positive
    if (! all(diff(ubreaks)>0))
      stop ("break points invalid: length of histogram cells must be greater than 0")
  }

  ## number of bins
  nbins <- length(ubreaks)-1

  if (nbins > 1e9)
    stop ("Argument 'breaks': too many bins (> 1e9).")

  ## --- pre-sample and type of distribution --------------------------------

  ## random sample of size n
  ## ("pre-sample" required to get some information about the distribution)
  X <- myrdist(n)

  ## estimate distribution type
  if (is.null(dtype)) {
    tmp <- X[1:100]
    if (isTRUE(all.equal(tmp,round(tmp))) && isTRUE(all(tmp<1e6))) {
      dtype <- "discr"
    } else {
      dtype <- "cont"
    }
  }

  ## check given data again
  if (dtype=="discr") {
    if(isTRUE(exactu)) {
      warning("Argument 'exactu' ignored for discrete distributions.")
      exactu <- FALSE
    }
    if (is.null(pdist)) {
      stop ("Argument 'pdist' required for discrete distribution.")
    }
    if (!is.null(qdist)) {
      warning("Argument 'qdist' ignored for discrete distributions.")
      qdist <- NULL
    }
  }
  
  ## --- 'qdist' given: compute break points in x-scale ---------------------

  if (! is.null(qdist)) {
    xbreaks <- qdist(ubreaks,...)
  }
  
  ## --- 'qdist' not given: recompute break points in u- and x-scale --------
  
  if (is.null(qdist)) {

    if (! isTRUE(exactu)) {
      ## it is faster to have break points in x-scale.
      ## if allowed we use the empirial quantiles of the first sample
      ## and recompute the break points in u-scale.

      ## compute empirial quantiles
      xbreaks <- quantile(X, probs=ubreaks, na.rm=TRUE)
      names(xbreaks) <- NULL
      xbreaks[1]       <- if (is.null(trunc)) -Inf else trunc[1]
      xbreaks[nbins+1] <- if (is.null(trunc))  Inf else trunc[2]

      ## adjust break points in u-scale
      ubreaks <- pdist(xbreaks,...)
      ubreaks[ubreaks<0] <- 0
      ubreaks[ubreaks>1] <- 1
      
    } else {
      ## break points in x-scale not available
      xbreaks <- rep(NA,length(ubreaks))
    }
  }
    
  ## --- check probabilities ------------------------------------------------

  ## we have problems when expected probabilities are too small.
  ## however, this may happens for discrete random variates.
  ## thus we collapse the corresponding bins.

  ## expected probabilities
  p0 <- diff(ubreaks)

  ## find bins which are "too small"
  too.small <- which(p0<min.bin.width)
  
  if (length(too.small)>0) {
    ## collapse bins
    if (dtype != "discr")
      warning("probability for some bins too small --> collapse bins.")

    ubreaks <- ubreaks[-too.small]
    ubreaks[length(ubreaks)] <- 1

    xbreaks.max <- xbreaks[length(xbreaks)]
    xbreaks <- xbreaks[-too.small]
    xbreaks[length(xbreaks)] <- xbreaks.max

    nbins <- length(ubreaks)-1
  }
  
  ## --- compute frequency tables -------------------------------------------

  ## table for storing frequencies
  count <- matrix(0,nrow=rep,ncol=nbins)
     
  ## loop for each row of table
  for (i in 1:rep) {

    ## random sample of size n
    ## (for i==0 we reuse the pre-sample)
    if (i>1) X <- myrdist(n)
    ## X must be of class "numeric"
    if (is.integer(X)) { X <- as.numeric(X) }

    ## get row
    if (! is.na(xbreaks[1])) {
      ## we can construct the histogram using the x-values
      count[i,] <- .Call("rvgt_bincount",X,xbreaks,PACKAGE="rvgtest")
    }
    else {
      ## otherwise we first have to transform the x-values into
      ## uniformly distributed u-values
      ## (slower but more robust for densities with poles)

      ## get frequency table (using function hist)
      u <- pdist(X,...)
      count[i,] <- .Call("rvgt_bincount",u,ubreaks,PACKAGE="rvgtest")
    }
  }
        
  ## --- prepare result -----------------------------------------------------

  ## return result as object of class "rvgt.ftable"
  ftable <- list(n=n,rep=rep,ubreaks=ubreaks,xbreaks=xbreaks,count=count,dtype=dtype)
  class(ftable) <- "rvgt.ftable"

  ## plot histogram
  if( isTRUE(plot) ) {
    plot(ftable)
  }

  ## return frequency table
  return(invisible(ftable))
}

## --------------------------------------------------------------------------

invqdist <- function (x, qdist, ...)
  ## ------------------------------------------------------------------------
  ## Compute distribution function by inverting quantile function.
  ## ------------------------------------------------------------------------
  ## x      : Point where CDF has to be evaluated
  ## qdist  : Quantile function of distribution
  ## ....   : Parameters of distribution
  ## ------------------------------------------------------------------------
{
  ql <- qdist(0,...); if (isTRUE(ql>=x)) return(0)
  qu <- qdist(1,...); if (isTRUE(qu<=x)) return(1)

  f <- function(z) { qdist(z,...) - x } 

  return (uniroot(f, interval=c(0,1), f.lower=(ql-x), f.upper=(qu-x),
                  tol=1e-100 * .Machine$double.eps)$root)
}

## --------------------------------------------------------------------------
