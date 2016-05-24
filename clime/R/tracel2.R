tracel2 <- function(Sigma, Omega) {
  return(sum(diag(  (Sigma%*%Omega -  diag(1,  dim(Omega)[1] ) )^2 )) )
}
