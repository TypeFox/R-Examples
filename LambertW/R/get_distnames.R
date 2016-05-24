#' @rdname distname-utils
#' 
#' @description 
#' \code{get_distnames} lists all currently implemented distributions \eqn{F_X}.
#' 
#' @return 
#' \code{get_distnames} returns a vector of strings in alphabetical order. 
#' It lists all supported distributions. 
#' Each string can be passed as the \code{distname} argument to several functions in this package.
#' @export

get_distnames <- function() {
  # keep them in alphabetical order; otherwise tests fail
  implemented.dists <- c("cauchy", 
                         "chisq", 
                         "exp",
                         "f",  # for F distribution use 'f' to avoid confusion with 'F' for FALSE in R
                         "gamma", 
                         "laplace", 
                         "normal", 
                         "t", 
                         "unif")
  return(implemented.dists)
}