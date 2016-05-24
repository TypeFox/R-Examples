raise.matrix <- 
function(x, power = 1)
{
   if (!inherits(x, "matrix") || nrow(x) != ncol(x))
      stop("'x' must be a square matrix!")
   if (length(power) != 1)
      stop("'power' must be a single real number!")
   eig <- eigen(x)
   val <- diag(eig$values ^ power)
   vec <- eig$vectors
   out <- vec %*% val %*% t(vec)
   return(out)
}
