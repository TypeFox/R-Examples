#' Enumerate the partitions of an integer
#' 
#' Determine all unrestricted partitions of an integer.  This
#' function is equivalent to the function parts in the partitions
#' package.
#' 
#' @param n an integer
#' @return a matrix whose rows are the n-tuples
#' @author Robin K. S. Hankin, from package partitions
#' @export
#' @examples
#' partitions(5)
#' str(partitions(5))
partitions <- function(n){
  stopifnot(all(sapply(as.list(n), is.wholenumber)))
  stopifnot(all(n >= 1))
  if(is.vector(n) && length(n) > 1){
    return(lapply(as.list(n), function(x){
      as.matrix( t(parts(x)) )
    }))
  }
  as.matrix( t(parts(n)) )
}