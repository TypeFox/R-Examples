toeplitz_nonsym <- function(vec){
	N <- length(vec)
	mat <- matrix(0,2*N-1,N)
	for(i in 1:N){
		mat[,i] <- c(rep(0,i-1),vec,rep(0,N-i))
	}

	return(mat)
}

scaling_function <- function(filter,J){

## Computes the scaling function and the wavelet function (for compactly supported wavelet) using the
## cascade algorithm on the grid of dyadic integer 2^{-J}
##
## 	INPUT	filter	Wavelet filter
##        	J   	Largest scale
##
## 	OUTPUT  phi 	Scaling function
##         	psi 	Wavelet function
##
##					based on the paper of Fay, Moulines, Roueff, Taqqu, 2009
##                                      Achard & Gannaz (2014)
##________________________________________________________________________________________________

##
##... compute first the scaling functions on the integers
##
N <- length(filter)
Htmp <- toeplitz_nonsym(filter)
H <- matrix(0,N,N)
for(irow in 1:N){
	H[irow,] <- Htmp[2*(irow-1)+1,]
} 
H <- sqrt(2)*H 
	res_eigen <- eigen(H)
	V <- res_eigen$vectors
	D <- diag(res_eigen$values)

j <- 1
while((abs(D[j,j]-1) >= 1.0e-5)&(j<N)){
    j <- j+1
}
phi <- V[,j]
phi <- phi/sum(phi)

##
##... cascade algorithm
##
a <- sqrt(2)*filter
for(j in 1:(J-1)){
    phi <- convolve(phi,rev(a),type='open')
    a <- as.vector(rbind(a,rep(0,N)))
} 
   
##
##... compute the final scaling / wavelet function (using the wavelet
##... scaling equation)
##
b <- sqrt(2)*rev(filter*(-1)^(0:(N-1)))	
b <- as.vector(rbind(b,array(0,dim=c(2^(J-1)-1,N))))
psi <- convolve(phi,rev(b),type='open')
phi <- convolve(phi,rev(a),type='open')

phi <- phi[1:(2^J*(N-1)+1)]
psi <- psi[1:(2^J*(N-1)+1)]

list(phi=phi,psi=psi)

}
