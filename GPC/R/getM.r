getM <- function(d,p){
  M <- choose(d+p,d)
  # M <- factorial(N+P)/factorial(N)/factorial(P)
	return(M)
}