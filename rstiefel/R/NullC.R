NullC <-
function(M)
{
  #modified from package "MASS" 
  #MASS version : Null(matrix(0,4,2))  returns a 4*2 matrix
  #this version : NullC(matrix(0,4,2)) returns diag(4)

  tmp <- qr(M)
  set <- if (tmp$rank == 0L)
      1L:nrow(M)
  else -(1L:tmp$rank)
  qr.Q(tmp, complete = TRUE)[, set, drop = FALSE]
}
