#' @title Simple Random Sample
#' 
#' @description Computes all possible samples from a given population using simple random sampling
#' 
#' @details If non-finite values are entered as part of the population, they are removed; and the returned simple random sample computed is based on the remaining finite values. 
#' 
#' @param popvalues are values of the population.  \code{NA}s and \code{Inf}s are allowed but will be removed from the population. 
#' @param n the sample size 
#' 
#' @return The function \code{srs()} returns a matrix containing the possible simple random samples of size \code{n} taken from a population of finite values \code{popvalues}.
#' 
#' @author Alan T. Arnholt <arnholtat@@appstate.edu> 
#' 
#' @seealso \code{\link{combn}}
#' @export
#'  
#' @examples
#' 
#' srs(popvalues = c(5, 8, 3, NA, Inf), n = 2)
#'
#' @keywords programming
#####################################################################################
# Updated 6/18/13 to work/exclude non-finite values--and to use combn() 
# instead of homegrown Combinations()
#####################################################################################
srs <- function(popvalues, n){
popvalues <- popvalues[is.finite(popvalues)]
N <- length(popvalues)
store <- t(combn(N, n))
matrix(popvalues[t(store)], nrow = nrow(store), byrow = TRUE)
}
