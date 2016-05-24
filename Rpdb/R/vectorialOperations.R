## Vectorial operations

dotProd <- function(U,V){
  if(!(is.vector(U) & length(U)==3 & is.numeric(U)))
    stop("'U' must be a numeric vector of length 3")
  if(!(is.vector(V) & length(V)==3 & is.numeric(V)))
    stop("'V' must be a numeric vector of length 3")
  sum(U*V)
}

vectNorm <- function(U){
  if(!(is.vector(U) & length(U)==3 & is.numeric(U)))
    stop("'U' must be a numeric vector of length 3")
  sqrt(dotProd(U,U))
}

rotVect <- function(U, n = 1){
  if(!(is.vector(U) & length(U)==3 & is.numeric(U)))
    stop("'U' must be a numeric vector of length 3")
  if(!(length(n)==1 & round(n)==n))
    stop("'n' must be an integer")
  if(n>0)
    U <- U[c((length(U)-n+1):length(U),1:(length(U)-n))]
  if(n<0)
    U <- U[c((abs(n)+1):length(U),1:abs(n))]
  return(U)
}

vectProd <- function(U,V){
  if(!(is.vector(U) & length(U)==3 & is.numeric(U)))
    stop("'U' must be a numeric vector of length 3")
  if(!(is.vector(V) & length(V)==3 & is.numeric(V)))
    stop("'V' must be a numeric vector of length 3")
  rotVect(U,-1) * rotVect(V, 1) - rotVect(U, 1) * rotVect(V,-1)
}
  

