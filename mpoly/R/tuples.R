#' Determine all n-tuples using the elements of a set.
#' 
#' Determine all n-tuples using the elements of a set.  This is
#' really just a simple wrapper for expand.grid, so it is not
#' optimized.
#' 
#' @param set a set
#' @param n length of each tuple
#' @param repeats if set contains duplicates, should the result?
#' @param list tuples as list?
#' @return a matrix whose rows are the n-tuples
#' @export
#' @examples
#' 
#' tuples(1:2, 5)
#' tuples(1:2, 5, list = TRUE)
#' 
#' apply(tuples(c("x","y","z"), 3), 1, paste, collapse = "")
#' 
#' # multinomial coefficients
#' r <- 2 # number of variables, e.g. x, y
#' n <- 2 # power, e.g. (x+y)^2
#' apply(burst(n,r), 1, function(v) factorial(n)/ prod(factorial(v))) # x, y, xy
#' mp("x + y")^n
#' 
#' r <- 2 # number of variables, e.g. x, y
#' n <- 3 # power, e.g. (x+y)^3
#' apply(burst(n,r), 1, function(v) factorial(n)/ prod(factorial(v))) # x, y, xy
#' mp("x + y")^n
#' 
#' r <- 3 # number of variables, e.g. x, y, z
#' n <- 2 # power, e.g. (x+y+z)^2
#' apply(burst(n,r), 1, function(v) factorial(n)/ prod(factorial(v))) # x, y, z, xy, xz, yz
#' mp("x + y + z")^n
#' 
#' 
tuples <- function(set, n = length(set), repeats = FALSE, list = FALSE){

  ## determine how big the output will be
  r <- length(set)	
  nCombos <- r^n
  
  ## expand.grid really does all the work in this function, so
  ## setup the interior part of expand.grid.  
  ## sets4call looks like "set2, set2, set2" with n = 3
  out <- do.call(expand.grid, replicate(n, set, simplify = FALSE))
  out <- unname(as.matrix(out))
  out <- out[,n:1]
  
  ## delete duplicates
  if(!repeats) out <- unique(out)

  ## do list
  if(list) out <- unname(split(out, rep(1:nCombos, n)))
  
  ## return
  out
}

