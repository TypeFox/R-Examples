compute_nj <- function(n,N){

## Computes the number of wavelet coefficients at each scale
##
## 	INPUT	n  	Sample size
##        	N   	Filter length
##
## 	OUTPUT  nj	Number of wavelet coefficients at each scale
##         	J	Number of scales
##
##					based on the paper of Fay, Moulines, Roueff, Taqqu, 2009
##                                      Achard & Gannaz (2014)
##________________________________________________________________________________________________

nj <- vector()

nj[1] <- n
J <- 2
while((nj[J-1]-N +1) >= 1){
	nj[J] <- ceiling((nj[J-1]-N+1)/2)
	J <- J+1
}
J <- J-2
nj <- nj[2:length(nj)]

list(nj=nj,J=J)
}
