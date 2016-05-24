#' Classic ALK
#' 
#' Returns an Age-Length Key calculated from a matrix with the count of
#' individuals per age- and length-class, as described by Fridriksson (1934).
#' 
#' @param x A matrix with \code{i} lines and \code{j} columns, where
#' \code{x[i, j]} is the count of individuals of length \code{i} and age
#' \code{j}.
#' 
#' @return A matrix with the probability of an individual of length \code{i}
#' having age \code{j}, i.e. \eqn{P(j|i)}.
#' 
#' @references Fridriksson, A (1934). On the calculation of age-distribution
#' within a stock of cod by means of relatively few age determinations as a key
#' to measurements on a large scale. \emph{Rapp. P.-V. CIEM}, \strong{86}, 1-5.
#' 
#' @examples
#' data(hom)
#' calc_ALK(hom$otoliths[[1]])
#' 
#' @export
calc_ALK <- function(x) {
  ni <- rowSums(x, na.rm = TRUE)
  ni[ni == 0] <- 1
  return(x / ni)
}

