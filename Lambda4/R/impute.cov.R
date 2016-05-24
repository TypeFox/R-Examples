#' Compute Covariance Matrix
#' 
#' @description Implements various missing data techniques and generates a covariance matrix.
#' 
#' @param x A data matrix
#' @param missing how to handle missing values.
#' 
#' @author Tyler Hunt \email{tyler@@psychoanalytix.com}
#' 
#' @export

impute.cov <- function(x, missing = c('complete', 'pairwise', 'mi'))
{
  p <- ncol(x)
  if (nrow(x) == p)
    x <- as.matrix(x)
  else{
    missing <- match.arg(missing)
    switch(missing
           , complete = cov(x, use = "complete")
           , pairwise = cov(x, use = "pairwise")
           , mi = {
             stopifnot(require(mice))
             cov(complete(mice(x, diagnostics = FALSE)))
           }
    )
  }
}
