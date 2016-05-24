my_check <-
function(y, X, perm)
{
  ## Internal function to check argument:s y, X, perm
  
  # y as numeric vector
  if (!is.vector(y) || mode(y) != "numeric")
    stop("Sorry, argument 'y' must be a numeric vector")
  # no misssing values in y
  if (any(is.na(y))) 
    stop("No missing data allowed in argument 'y' ")	
  # binary values (0, 1) in y
  if (!all(y %in% c(0, 1)))
    stop("Argument 'y' must contain only 0 and 1")
  # X as matrix or data frame
  if(!is.matrix(X) & !is.data.frame(X))
    stop("Argument 'X' must be a matrix or data.frame")    
  # compatibility between X and y
  if (nrow(X) != length(y)) 
    stop("'X' and 'y' have different lengths")
  # force X as matrix
  if (!is.matrix(X)) X = as.matrix(X)
  # permutations
  if (mode(perm) != "numeric" || length(perm) != 1
      || perm < 0 || (perm %% 1) !=0) 
  {
    warning("argument 'perm' incorrectly defined. Value perm=100 is used")
    perm = 100
  }
  # results
  list(y=y, X=X, perm=perm)
}