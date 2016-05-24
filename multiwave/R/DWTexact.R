DWTexact <- function(x,filter){

## Computes the scaling function and the wavelet function (for compactly supported wavelet) using the
## cascade algorithm on the grid of dyadic integer $2^{-\text{J}}$ 
##
## 	INPUT	x		Vector of signal values to 
##		filter  	Filter of the wavelet transform
##
## 	OUTPUT  dwt		Vector of wavelet coeffcients
##         	indmaxband 	Vector containing the maximal index of wavelet coefficients at each scale
##		Jmax		Maximal scale
##
##					based on the paper of Fay, Moulines, Roueff, Taqqu, 2009
##                                      Achard & Gannaz (2014)
##________________________________________________________________________________________________

n <- length(x)
N <- length(filter)
h<-filter
g<-rev(filter*(-1)^(0:(N-1)))

tmp <- compute_nj(n,N)
nscale <- tmp$nj
Jmax <- tmp$J
indmaxband <- cumsum(nscale)

# Pyramidal algorithm
dwt <- rep(0,indmaxband[(length(indmaxband))])

aold <- x
indmin <- 1
for(iscale in 1:Jmax){
	acur <- rep(0,nscale[iscale])
 	d <- rep(0,nscale[iscale])
	offset <- 1
	for(p in 1:nscale[iscale]){
		acur[p] <- sum(h*aold[offset:(offset+N-1)])
		dwt[indmin+p-1] <- sum(g*aold[offset:(offset+N-1)])
		offset <- offset+2
  	}
  	aold <- acur
  	indmin <- indmaxband[iscale]+1
}


list(dwt=dwt,indmaxband=indmaxband,Jmax=Jmax)


}
