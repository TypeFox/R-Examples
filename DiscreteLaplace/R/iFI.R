iFI <- function(p,q)
{
M <- matrix(p * q * (1 - p) * (1 - q)/(1 + p * q), 2, 2)
M[1, 1] <- M[1, 1] * (1 - p) * (1 - p * q^2)/q/(1 - q)^2
M[2, 2] <- M[2, 2] * (1 - q) * (1 - q * p^2)/p/(1 - p)^2
return(M)
}
