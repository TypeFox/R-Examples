#needed for nonmetric sphere primal

vecAsDist <- function(x)
{
  n <- (1+sqrt(1+8*length(x)))/2
  e <- matrix(0,n,n); k<-0
  for (i in 1:(n-1)) {
    l <- n-i
    ll <- 1:l
	  e[i+ll,i] <- x[k+ll]
   k<-k+l
	}
return(as.dist(e))
}