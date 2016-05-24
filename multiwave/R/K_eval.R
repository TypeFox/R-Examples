K_eval <- function(psi_hat,u,d){
## computes the matrix whith elements
## K(dl,dm)=int( u^(dl+dm) |psi_hat(u)|^2 du )
## 
## 	INPUT	 d: kx1 parameter value
##               psi_hat: the Fourier transform of the wavelet mother at
##               values u
##               u: the grid for the approximation of the integral
##
##                                           Achard & Gannaz (2014)
##______________________________________________________________________

l <- length(psi_hat)
k <- length(d)

K_fct <- matrix(0,k,l)

for(j in 1:k){
    K_fct[j,] <- (abs(u+(u==0))^(-d[j])-(u==0))*psi_hat
}
K <- K_fct%*%t(Conj(K_fct))*(max(u)-min(u))/l

K <- Re(K) ### K is real by construction

return(K)
}
