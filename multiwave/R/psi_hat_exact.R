psi_hat_exact<-function(filter,J){

## Computes the discrete Fourier transform of the wavelet associated to 
## the given filter
##
## The wavelet is given by the function makescalingfunction of Fay et al
## (2009).
##
## The length of the Fourier transform is equal to the length of the grid
## where the wavelet is evaluated
##
## 	INPUT	filter  Quadratic mirror filter of the wavelet
##              J	Number of scales where the wavelet is evaluated
##
##	OUTPUT	fct	Values of the discrete Fourier transform of the wavelet
##		grid	Frequencies where the Fourier transform is evaluated
##
##                                           		Achard & Gannaz (2014)
##________________________________________________________________________________

Jmax <- J-3
res <- scaling_function(filter,J)
phi <- res$phi
psi <- res$psi
L <- length(psi)
psi <- psi/sqrt(2^J)

Nfft <- 2^round(log(L)/log(2))
Te <- (length(filter)-1)*2^Jmax ## support-size
psi <- c(psi,rep(0,(Nfft-L)))
psih <- fft(psi)/sqrt(pi*Te)

grid <- pi*Te*c(seq((1/Nfft),1/2,(1/Nfft)),seq((-1/2),(-1/Nfft),(1/Nfft)))


list(psih=psih,grid=grid)

}

