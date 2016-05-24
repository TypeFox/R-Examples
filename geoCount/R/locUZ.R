
####################################
#### Uniformize locations
####################################
unifLoc <- function(loc, length=1){
    lmax <- max(loc); lmin <- min(loc)
    res <- (loc - lmin)/(lmax - lmin) 
    res*length
  }
####################################
#### locations to Distance matrix
####################################
loc2U <- function(loc){
  if(is.data.frame(loc)) loc <- as.matrix(loc)
  if(!is.matrix(loc)) {
    stop("loc needs to be a matrix or data.frame!")
    return()
  }
  .Call( "loc2Ucpp", loc, PACKAGE = "geoCount" )
}
# R version
loc2U_R <- function(loc){
  if(is.data.frame(loc)) loc <- as.matrix(loc)
  if(!is.matrix(loc)) {
    stop("loc needs to be a matrix or data.frame!")
    return()
  }
n <- nrow(loc)
dst <- function(t){sqrt((t[1]-t[3])^2+(t[2]-t[4])^2)}
U <- matrix(0, n, n)
for(i in 1:n){
  U[,i] <- apply(cbind(matrix(rep(loc[i,],each=n),n,),loc),1,dst)
}
U
}
####################################
#### locations to Distance matrix
####################################
locUloc <- function(loc, locp){
  if(is.data.frame(loc)) loc <- as.matrix(loc)
  if(is.data.frame(locp)) locp <- as.matrix(locp)
  if(is.vector(loc)) loc <- matrix(loc,,2)
  if(is.vector(locp)) locp <- matrix(locp,,2)
  if( (!is.matrix(loc)) || (!is.matrix(locp)) ){
    stop("Both loc and locp needs to be a matrix or data.frame!")
    return()
  } 
  .Call( "locUloccpp", loc, locp, PACKAGE = "geoCount" )
}
# R version
locUloc_R <- function(loc, locp){
  if(is.data.frame(loc)) loc <- as.matrix(loc)
  if(is.data.frame(locp)) locp <- as.matrix(locp)
  if(is.vector(loc)) loc <- matrix(loc,,2)
  if(is.vector(locp)) locp <- matrix(locp,,2)
  if( (!is.matrix(loc)) || (!is.matrix(locp)) ){
    stop("Both loc and locp needs to be a matrix or data.frame!")
    return()
  }
  n <- nrow(loc); np <- nrow(locp)
  dst <- function(t){sqrt((t[1]-t[3])^2+(t[2]-t[4])^2)}
  res <- matrix(0, np,n )
  for(i in 1:n){    
    res[,i] <- apply(cbind(matrix(rep(loc[i,], each=np),np,), locp),1,dst)
  }
  res
}
####################################
#### Correlation functions
####################################
rhoPowerExp <- function(u, a, k) { exp(-(u/a)^k ) }
rhoSph <- function(u, a, k=NULL) { 
  ifelse(u < a, 1 - 1.5*(u/a) + 0.5*(u/a)^3, 0) 
}
rhoMatern <- function(u, a, k){
  phi <- a; kappa <- k
  res <- ifelse(u > 0,   ( 2^(1-kappa) / gamma(kappa) ) * (u/phi)^kappa * 
besselK(x = u/phi, nu=kappa, expon.scaled = FALSE) , 1)
  res
}
#########################################
#### Distance matrix to Covariance matrix
#########################################
U2Z <- function(U, cov.par, rho.family = "rhoPowerExp"){
  s <- cov.par[1]; a <- cov.par[2]; k <- cov.par[3]
  if(rho.family=="rhoPowerExp"){
      Z <- s^2* rhoPowerExp(U, a, k)
    } else if(rho.family=="rhoMatern"){
        Z <- s^2* rhoMatern(U, a, k)
      } else if(rho.family=="rhoSph"){
          Z <- s^2* rhoSph(U, a, k)
      } else {
          cat("Notice: rho.family=", rho.family, " is not appropriate! rho.family=rhoPowerExp will be used.\n", sep="")
          Z <- s^2* rhoPowerExp(U, a, k)
        }
  Z
}
####################################
#### END
####################################