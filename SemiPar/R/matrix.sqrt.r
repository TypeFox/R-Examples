########## S-function: matrix.sqrt ##########

# Computes the positive square-root of a square
# positive definite matrix.

# Last changed: 05/06/97

matrix.sqrt <- function(A)
{
   sva <- svd(A)
   if (min(sva$d)>=0)
      Asqrt <- t(sva$v %*% (t(sva$u) * sqrt(sva$d)))
   else
      stop("Matrix square root is not defined")
   return(Asqrt)
}

######### End of S-function matrix.sqrt ########
