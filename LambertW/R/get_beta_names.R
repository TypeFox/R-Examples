#' @rdname beta-utils
#' @description
#' \code{get_beta_names} returns (typical) names for each component of
#' \eqn{\boldsymbol \beta}.
#' 
#' Depending on the distribution 
#' \eqn{\boldsymbol \beta} has different length and names: e.g., 
#' for a \code{"normal"} distribution \code{beta} is of length 
#' \eqn{2} (\code{"mu"}, \code{"sigma"}); for an \code{"exp"}onential 
#' distribution \code{beta} is a scalar (rate \code{"lambda"}).
#' 
#' @return
#' \code{get_beta_names} returns a vector of characters.
#' @export
#' 
get_beta_names <- function(distname) {
  check_distname(distname)
  
  switch(distname,
         cauchy = {beta.names <- c("location", "scale")}, # cauchy does not have finite mean/variance          
         chisq = {beta.names <- "df"},
         exp = {beta.names <- c("lambda")},
         f = {beta.names <- c("df1", "df2")},
         gamma = {beta.names <- c("shape", "scale")},
         laplace = {beta.names <- c("location", "scale")},
         normal = {beta.names <- c("mu", "sigma")},
         t = {beta.names <- c("location", "scale", "df")},
         unif = {beta.names <- c("min", "max")},
         user = {beta.names <- NULL})  # user defined distributions return NULL
  return(beta.names)
} 
