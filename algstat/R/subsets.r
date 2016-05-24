#' Compute Subsets
#'
#' Compute the subsets of a given set.
#'
#' Note that this algorithm is run in R: it is therefore not intended to be the most efficient algorithm for computins subsets.
#' 
#' @param set the original set
#' @param sizes desired size(s) of subsets
#' @param include_null should the empty vector be included?
#' @return a list of subsets as vectors
#' @export subsets
#' @seealso \code{\link{combn}}
#' @examples
#' 
#' 
#' subsets(1:3)
#' subsets(1:3, size = 2)
#' subsets(1:3, include_null = TRUE)
#' 
#' subsets(c('a','b','c','d'))
#' subsets(c('a','b','c','d'), include_null = TRUE)
#' 
#'
#' 
#' 
#'
#' 
#' 
subsets <- function(set, sizes = 1:length(set), include_null = FALSE){
  if((length(set) == 1) && (is.numeric(set) || is.integer(set) )){
    set <- 1:set
  }
  subsetsBySize <- lapply(sizes, function(n){  
    combn(length(set), n, function(x){
  	  set[x]
    }, simplify = FALSE)
  })
  out <- unlist(subsetsBySize, recursive = FALSE)
  if(include_null) out <- unlist(list(list(set[0]), out), recursive = FALSE)
  out
}