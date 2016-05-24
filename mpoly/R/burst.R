#' Enumerate integer r-vectors summing to n
#' 
#' Determine all r-vectors with nonnegative integer entries summing 
#' to n.  Note that this is intended to be optimized.
#' 
#' @param n integer to sum to
#' @param r number of components
#' @return a matrix whose rows are the n-tuples
#' @export
#' @examples
#' burst(4)
#' 
#' burst(4, 4)
#' burst(4, 3)
#' burst(4, 2)
#' 
#' rowSums(burst(4))
#' rowSums(burst(4, 3))
#' rowSums(burst(4, 2))
#' 
#' 
#' burst(10, 4) # all possible 2x2 contingency tables with n=10
#' burst(10, 4) / 10 # all possible empirical relative frequencies
#' 
burst <- function(n, r = n){
  stopifnot(is.wholenumber(n))
  stopifnot(n > 0)
  stopifnot(is.wholenumber(r))
  stopifnot(r > 0)
  
  if(r == 1) return(n)
   
  ## compute all partitions of the number (order does not matter)
  parts <- partitions(n)

  ## convert to list, remove 0's
  partsWOzeros <- apply(parts, 1, function(row) row[row != 0] )
  
  ## select all those with length less than or equal to r
  isAtLongest <- function(howLong){
    function(v) length(v) <= howLong
  }
  #isAtLongest(4)(1:4)
  #isAtLongest(4)(1:5)
  partsWOzeros <- Filter(isAtLongest(r), partsWOzeros)
  
  ## refill with 0's
  rvectors <- lapply(partsWOzeros, function(v){
    if(length(v) < r) return(c(v, rep(0, r - length(v))))
    v
  })
    
  ## compute permutations of each
  permsOfEach <- lapply(rvectors, permutations)
  
  ## make out matrix
  m <- sum(vapply(permsOfEach, nrow, integer(1)))
  n <- r
  out <- integer(length = m*n)
  dim(out) <- c(m, n)
  
  ## populate out matrix
  row <- 1
  for(i in seq_along(permsOfEach)){
    permsInPartition <- nrow(permsOfEach[[i]])
    out[row:(row+permsInPartition-1),] <- permsOfEach[[i]]
    row <- row + permsInPartition
  }
  
  ## return
  out
}









