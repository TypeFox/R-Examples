#' @rdname distname-utils
#' 
#' @description
#' \code{get_distname_family} determines whether a distribution is a 
#' location, scale, or location-scale family. 
#' It also returns whether the distribution is supported on non-negative
#' values only.
#' 
#' @return 
#' \code{get_distname_family} returns a list with
#' \item{location}{ logical; if \code{TRUE}, the distribution is a location family,}
#' \item{scale}{ logical; if \code{TRUE}, the distribution is a scale family.}
#' \item{is.non.negative}{ logical; if \code{TRUE}, the distribution has support only for
#' the non-negative reals (this is usually the case when \code{location = FALSE}
#' and \code{scale = TRUE})}

#' @export
#' @examples
#' 
#' get_distname_family("normal")


get_distname_family <- function(distname) {
  
  check_distname(distname)
  out <- list(location = FALSE,
              scale = FALSE,
              is.non.negative = FALSE)
  if (distname %in% c("cauchy", "normal", "t", "unif")) {
    out$location <- TRUE
    out$scale <- TRUE
  } else if (distname %in% c("exp", "chisq", "gamma", "F", "f")) {
    out$scale <- TRUE
    out$is.non.negative <- TRUE
  }
  return(out)
}