#' @encoding UTF-8
#' @title The Missing R Parallel Sum
#'
#' @description Provides parallel sum like \code{pmin} and \code{pmax} from the base package. The function \code{sum} simply does not help when the objective is to obtain a vector with parallel sum rather than a scalar value.
#'
#' @param \dots One or more unit objects
#' @param na.rm A logical value \code{TRUE} or \code{FALSE}, the default
#'
#' @return A vector containing the parallel sum.
#'
#' @author Daniel Marcelino, \email{dmarcelino@@live.com}
#'
#' @keywords Misc
#'
#' @examples
#' if (interactive()) {
#' n <- 20;
#' mydat <- data.frame(PT = rnorm(n, mean = .30),
#' PSDB = rnorm(n, mean = .25), PSB = rnorm(n, mean = .15));
#' head(mydat);
#' transform(mydat, DK = psum(PT, PSDB, PSB - 1));
#' }
#'
#' @export
`psum` <-
  function(..., na.rm=FALSE) {
    x <- list(...)
    rowSums(matrix(unlist(x), ncol=length(x)), na.rm=na.rm)
  }
NULL
