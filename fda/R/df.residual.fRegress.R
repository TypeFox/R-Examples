df.residual.fRegress <- function(object, ...){
  dfm <-   object$df
  if(is.null(dfm))
    stop("'object' does not have a 'df' component.")
#
  nobs <- length(object$wt)
  dfr <- nobs-dfm
  attr(dfr, 'nobs') <- nobs
  attr(dfr, 'df.model') <- dfm
  dfr
}
