##Create two vectors of equal length
v1 <- rnorm(10)
v2 <- rnorm(10)

##compute the inner product between the vectors
dotProd(v1,v2)
##or
sum(v1*v2)

##compute the square 2-norm of v1
norm2(v1)
##or
dotProd(v1,v1)
##or
sum(v1*v1)

##If the vectors are of unequal length the longer vector
##gets truncated (with a warning). 
dotProd(v1,c(v2,2))
\dontshow{
  if( abs(dotProd(v1,v2)-sum(v1*v2)) > 1e-10 ){
    stop("dotProd: Results not equal")
  }
  if( (abs(norm2(v1) - dotProd(v1,v1)) > 1e-10) ||
     (abs(norm2(v1) - sum(v1*v1)) > 1e-10) ){
    stop("norm2: Results not equal")
  }
}
