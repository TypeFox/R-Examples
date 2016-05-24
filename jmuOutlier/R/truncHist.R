truncHist <-
function(x, xmin=NULL, xmax=NULL, trim=0.025, ...) {
  # Produces a truncated histogram.
  # This function may be useful if data contain some extreme outliers.
  # 'x': Vector of numerical observations.
  # 'xmin': Minimum numerical value to shown in graph. 
  # 'xmax': Maximum numerical value to shown in graph. 
  # 'trim': The fraction (0 to 0.5) of observations to be trimmed from
  #    each end of \code{x} before the histogram is constructed.
  #    The argument \code{trim} is used only when \code{xmin} and \code{xmax} are \code{NULL}.
  # ...: Optional arguments to \code{\link[graphics]{hist}}.
  # example:  truncHist( rnorm(1000) )
  # example:  truncHist( rcauchy(1000) )
  if (!is.numeric(x))  stop("'x' must be numeric.")
  if (length(xmin)>1 | (!is.numeric(xmin) & !is.null(xmin)))   stop("'xmin' must be scalar numeric or NULL.")
  if (length(xmax)>1 | (!is.numeric(xmax) & !is.null(xmax)))   stop("'xmax' must be scalar numeric or NULL.")
  if (length(trim)!=1 | !is.numeric(trim))   stop("'trim' must be scalar numeric between 0 and 0.5.")
  if (trim<0 | trim>0.5)   stop("'trim' must be scalar numeric between 0 and 0.5.")
  if ( (is.null(xmin)& !is.null(xmax)) || (!is.null(xmin)&is.null(xmax)))
     stop("'xmin' and 'xmax' must both be NULL or must both be nonNULL.")
  if (is.null(xmin)&is.null(xmax)) 
     y <- sort(x)[(length(x)*trim): (length(x)*(1-trim))]
  if (!is.null(xmin)&!is.null(xmax))  {
     if (xmin >= xmax) stop("'xmin' cannot be as large as or larger than 'xmax'.")
     y <- x[ x>xmin & x<xmax ]  }
  hist(y, ...)
}
