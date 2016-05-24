#' hist for ff vectors
#'
#' Currently this is a simple version of \code{\link{hist}} functionality.
#' @method hist ff 
#' @param x \code{ff} vector of values for which the histogram is desired
#' @param breaks a single numer given the number of cells for the histogram
#' @param plot logical. If \code{TRUE} (default), a histogram is plotted, otherwise a list of breaks and counts is returned
#' @param ... further arguments supplied to plot.
#' @return histogram object
#' @export
#' @export hist.ff
hist.ff <- function(x, breaks=min(100, length(x)), plot=TRUE, ...){
  xname <- deparse(substitute(x))
  
  # TODO improve breaks...
  if (length(breaks) <= 1){
    rng <- range(x)    
    breaks <- seq(from=rng[1], to=rng[2], length.out=breaks+1)
  }
  nB <- length(breaks)
  n <- length(x)
  
  equidist <- TRUE
  counts <- integer(nB-1)
  
  for (i in chunk(x, ...)){
    fi <- findInterval(x[i], breaks, rightmost.closed=TRUE)
    counts <- counts + tabulate(fi, nbins=nB-1)
  }
  
  dens <- counts/(n * diff(breaks))
  mids <- 0.5 * (breaks[-1L] + breaks[-nB])
  
  r <- structure( list( breaks = breaks
                      , counts = counts
                      , intensities = dens
                      , density = dens
                      , mids = mids
                      , xname = xname
                      , equidist = equidist
                      )
                , class = "histogram"
                )
  if (plot){
    plot(r, ...)
    invisible(r)
  } else {
    r
  }
}

#  x <- ff(rnorm(10000000))
#  hist.ff(x)
