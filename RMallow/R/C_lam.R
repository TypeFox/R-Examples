#' Calculate the normalizing coefficient for Mallow's model in a sequence
#' space.
#' 
#' Calculate the normalizing coefficient, as a function of the lambda
#' parameter, and the size of the sequence space.
#' 
#' 
#' @param lambda Spread parameter for Mallows' model.
#' @param dists Vector of all distances from each sequence to 1:N
#' @param dists.table Table version of "dists" above.
#' @return Normalizing coefficient of Mallows' model in N! space with lambda = lambda.
#' @author Erik Gregory
#' @keywords normalize

C_lam <-
function(lambda, dists = NULL, dists.table = NULL) {
  if (is.null(dists.table)) {
    dists.table <- table(dists)
  }
  ps <- dists.table*(exp(-lambda*as.numeric(names(dists.table))))
  C <- 1/sum(ps)
  return(C)
}
