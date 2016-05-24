BoxplotTrue <- function(x, truth=NULL, vnames=NULL, center=FALSE,
                         se.truth = NULL, color="white", ...) {
  ## Compares the marginal distributions of the columns of a matrix
  ## against a vector of true values.  Useful for validating an MCMC
  ## algorithm against simulated data where you know the true
  ## parameter values.
  ##
  ## Args:
  ##   x: A matrix of data to plot.  Rows are observations/MCMC
  ##     iterations.  Columns are variables.  Each column gets its own
  ##     box.
  ##   truth: A vector of "true" numeric values to compare against the
  ##     boxplots of x.
  ##   vnames:  A character vector to use as the boxplot labels.
  ##   center: Logical.  If TRUE then the value of 'truth' is
  ##     subtracted from each column of x before plotting.
  ##   se.truth: If the true standard errors from the simulation are
  ##     known then setting se.truth will add reference lines at +/- 2
  ##     standard errors.
  ##   color:  A vector of colors to use for the different boxes.
  ##   ...:  Extra arguments passed to boxplot.
  ##
  if(!is.null(truth)){
    if(length(truth)!=ncol(x))
      stop("truth and x don't conform in boxplot.true")
    if(!is.null(se.truth)){
      if(length(se.truth)!=length(truth)){
        stop("truth and se.truth don't conform in boxplot.true")
      }
    }
  }

  if(!is.null(se.truth) && is.null(truth)){
    stop("need to supply argument 'truth' to plot se.truth in boxplot.true")
  }

  if(center && !is.null(truth)){
    x <- x - matrix(rep(truth, nrow(x)), ncol=length(truth), byrow=TRUE)
    truth <- truth-truth
  }
  nx <- ncol(x)
  if(is.null(vnames)) vnames <- paste(1:nx)

  if(length(color) < nx) color <- rep(color, length.out=nx)
  boxplot(split(x, col(x)), names=vnames, col=color, ...)

  if(!is.null(truth)){
    AddSegments(1:ncol(x), truth, lwd=3)
    if(!is.null(se.truth)){
      AddSegments(1:ncol(x), truth+2*se.truth, lwd=1, lty=2)
      AddSegments(1:ncol(x), truth-2*se.truth, lwd=1, lty=2)
    }
  }
}
