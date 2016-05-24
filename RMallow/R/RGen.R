#' Initialize sequence modes for the clustering process.
#'
#' @param G number of cluster centers, including the hypothesis if provided
#' @param hyp a single sequence of length \code{abils} to initialize one of the cluster centers
#' @param abils number of items being ranked
#' @return A list of G cluster centers, each of length abils
#' @author Erik Gregory
#' @examples
#' Rgen(3, 1:5, 5)
Rgen <- function(G, hyp = NULL, abils) {
  R <- list()
  if (!is.null(hyp)) {
    R[[1]] <- hyp
  }
  else {
    R[[1]] <- sample(abils)
  }
  if (G > 1) {
    for (i in 2:G) {
      R[[i]] <- sample(abils)
    }
  }
  return(R)
}