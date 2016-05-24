# Transforms each column of a matrix using the functions in "trans"
.transMatrix <- function(mat, trans = NULL)
{
  parNames <- colnames(mat)
  
  # Transform the columns whose name is in "trans"
  if( !is.null(trans) )
  {
    trans <- unlist(trans)
    if( !all(sapply(trans, function(input) is.function(get(input)))) ) stop("All elements of \"trans\" should be names of functions")
    for(nam in names(trans))
    {
      index <- which( parNames == nam ) 
      mat[ , index] <- get(trans[nam])(mat[ , index])
    }
  }
  
  return( mat )
}