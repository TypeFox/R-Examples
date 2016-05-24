#' Element-wise arithmetic with vectors of multivariate polynomials.
#' 
#' Element-wise arithmetic with vectors of multivariate polynomials.
#' 
#' @param e1 an object of class mpolyList
#' @param e2 an object of class mpolyList
#' @return An object of class mpolyList.
#' @name mpolyListArithmetic
#' @examples
#' 
#' ( ms1 <- mp( c('x + 1', 'x + 2') ) )
#' ( ms2 <- mp( c('x + 1', 'y + 2') ) )
#' ms1 + ms2
#' ms1 - ms2
#' ms1 * ms2
#' 
NULL







#' @rdname mpolyListArithmetic
#' @export
`+.mpolyList` <- function(e1, e2){
  
  ## argument check
  if(is.numeric(e1) && length(e1) == 1) e1 <- mpoly(list(c(coef = e1)))
  if(is.mpoly(e1)) e1 <- mpolyList(e1)
  
  if(is.numeric(e2) && length(e2) == 1) e2 <- mpoly(list(c(coef = e2)))
  if(is.mpoly(e2)) e2 <- mpolyList(e2)
  
  if(!is.mpolyList(e1) || !is.mpolyList(e2)){
    stop('e1 and e2 must be of class mpolyList.', call. = FALSE)
  }
  
  if(length(e1) != length(e2)){
    stop('e1 and e2 must have equal length.', call. = FALSE)
  }
  
  
  ## determine length, flatten, and make indices on which to add
  n <- length(e1)
  
  flatList <- unlist(list(
    unclass(e1),
    unclass(e2)
  ), recursive = FALSE)
  
  ndcs2add <- split(cbind(1:n, (n+1):(2*n)), 1:n)
  
  
  ## sum
  out <- lapply(ndcs2add, function(v){
    Reduce('+', flatList[v])
  })
  out <- unname(out)
  
  
  ## caste and return
  class(out) <- 'mpolyList'
  out
}	












#' @rdname mpolyListArithmetic
#' @export
`-.mpolyList` <- function(e1, e2){
  ## change coefficient signs in e2
  e2 <- lapply(e2, unclass)
  e2 <- lapply(e2, function(l){
    lapply(l, function(v){
      v['coef'] <- -v['coef']
      v
    })	
  })
  class(e2) <- 'mpolyList'
    
  e1 + e2
}











#' @rdname mpolyListArithmetic
#' @export
`*.mpolyList` <- function(e1, e2){
  
  ## argument check
  
  if(is.numeric(e1) && length(e1) == 1) e1 <- mpoly(list(c(coef = e1)))  
  if(is.mpoly(e1)) e1 <- mpolyList(e1)
  
  if(is.numeric(e2) && length(e2) == 1) e2 <- mpoly(list(c(coef = e2)))  
  if(is.mpoly(e2)) e2 <- mpolyList(e2)
  
  stopifnot(is.mpolyList(e1))
  stopifnot(is.mpolyList(e2))
  
  if(length(e1) != length(e2)) stop('e1 and e2 must have equal length.', call. = FALSE)
  
  
  
  ## determine length, flatten, and make indices on which to multiply
  n <- length(e1)
  
  flatList <- unlist(list(
    unclass(e1),
    unclass(e2)
  ), recursive = FALSE)
  
  ndcs2add <- split(cbind(1:n, (n+1):(2*n)), 1:n)
  
  
  
  
  ## multiply
  out <- lapply(ndcs2add, function(v) Reduce('*', flatList[v]))
  out <- unname(out)
  
  
  
  ## caste and return
  class(out) <- 'mpolyList'
  out
}	







