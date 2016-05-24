########################### Random values from dirichlet distribution ###########################################
rdirich <- function(n,alphavec) {
 k <- length(alphavec)
 p <- matrix(rgamma(n*k,alphavec),n,k,byrow=T)
 sm <- matrix(apply(p,1,sum),n,k)
 return(p/sm)
}  
