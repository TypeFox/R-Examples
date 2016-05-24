#' Signed-Rank Estimate of Location (Intercept)
#' 
#' Returns the signed-rank estimate of intercept with is equivalent to the
#' Hodges-Lehmann estimate of the residuals.
#' 
#' 
#' @param x numeric vector
#' @return Returns the median of the Walsh averages.
#' @author John Kloke \email{kloke@@biostat.wisc.edu}
#' @seealso \code{\link{walsh}}
#' @references Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust
#' Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
#' 
#' Hollander, M. and Wolfe, D.A. (1999), \emph{Nonparametric Statistical
#' Methods}, New York: Wiley.
#' @examples
#' 
#' 
#' ## The function is currently defined as
#' function (x) 
#' median(walsh(x))
#' 
#' @export signedrank
signedrank <- function (x) median(walsh(x))
