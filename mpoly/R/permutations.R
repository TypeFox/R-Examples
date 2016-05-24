#' Determine all permutations of a set.
#' 
#' An implementation of the Steinhaus-Johnson-Trotter permutation
#' algorithm.
#' 
#' @param set a set
#' @return a matrix whose rows are the permutations of set
#' @export
#' @examples
#' permutations(1:3)
#' permutations(c('first','second','third'))
#' permutations(c(1,1,3))
#' apply(permutations(letters[1:6]), 1, paste, collapse = '')
permutations <- function(set){
  r <- length(set)
  if(r == 1 && is.numeric(set)) return(permutations(1:set))
    
  row2diag <- function(row, direction){
  	np1 <- length(row) + 1
  	mat <- matrix(nrow = np1, ncol = np1)
    if(direction == -1) for(k in 1:np1) mat[k,] <- insert(np1, np1-k+1, row)
    if(direction == +1) for(k in 1:np1) mat[k,] <- insert(np1, k, row)

    mat
  } 
  # row2diag(matrix(c(1,2,2,1),2,2)[,1], -1)
  # row2diag(matrix(c(1,2,2,1),2,2)[,1], +1)  
  
  stepUp <- function(mat){ # an r x # matrix
    c <- ncol(mat)
    m <- NULL
    for(k in 1:nrow(mat)) m <- rbind(m, row2diag(mat[k,],(-1)^k))
    m    
  }
  # stepUp(matrix(1))  
  # stepUp(matrix(c(1,2,2,1),2,2))
  
  # iterate stepUp
  out <- matrix(1)
  for(k in 1:(r-1)) out <- stepUp(out)
  
  # substitute set values
  for(k in 1:r) out[out==k] <- paste('_', set[k], sep = '')
  
  out <- gsub('_', '', out) # clear PH
  if(is.numeric(set)){
  	d <- dim(out)
    out <- as.numeric(out)
    dim(out) <- d
  }
  
  # reorder outcomes
  out <- out[do.call("order", split(out, col(out))),]
  
  # return
  unique(out)
}







#' Insert an element into a vector.
#'
#' Insert an element into a vector.
#' 
#' @param elem element to insert
#' @param slot location of insert
#' @param v vector to insert into
#' @return vector with element inserted
#' @export
#' @examples
#' insert(2, 1, 1)
#' insert(2, 2, 1) 
#' insert('x', 5, letters) 
insert <- function(elem, slot, v){
  n <- length(v)
  if(slot == 1) return( c(elem, v) )
  if(slot == n+1) return( c(v, elem) )  	
  c(v[1:(slot-1)], elem, v[slot:n])    
}