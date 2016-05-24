##--------------------------------------------------------------------------
## Author : Alexandre Brouste
## Project: Yuima
##--------------------------------------------------------------------------

##--------------------------------------------------------------------------
## Input                  mesh : mesh grid where the fBm is evaluated
##                        H    : self-similarity parameter 
##                        dim  : valued in R^dim
##      		 
##
## Output       :         simulation of a standard fractional Brownian noise
##                        for the mesh grid evaluation by Choleki s 
##			  decomposition of the covariance matrix of the fGn.
##--------------------------------------------------------------------------

##--------------------------------------------------------------------------
## Complexity O(N^3) via method chol.... mais mieux pour dim different
## Taille memoire N^2
## -------------------------------------------------------------------------


CholeskyfGn <- function( mesh , H , dim ){

	
H2 <- 2*H
mesh<-mesh[[1]]	
N<-length(mesh)-2 # N+1 is the size of the fGn sample to be simulated
fGn<-matrix(0,dim,N+1)

	
	
matcov <- matrix(0,N+1,N+1) # Covariance matrix of the fGn

		
for (i in (1:(N+1))) {
	j <- i:(N+1)
	matcov[i, j]<- 0.5 * (abs(mesh[i+1]-mesh[j])^H2 + abs(mesh[i]-mesh[j+1])^H2 - abs(mesh[i] - mesh[j])^H2-abs(mesh[i+1] - mesh[j+1])^H2)
	matcov[j, i]<- matcov[i, j]
}
	
L <- chol(matcov)
	
for (k in 1:dim){
	Z <- rnorm(N+1)
	fGn[k,] <- t(L) %*% Z
}	
		
return(fGn)

}