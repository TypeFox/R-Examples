mww_cov_eval<-function(d,x,filter,LU){

## Computes the multivariate wavelet Whittle estimator 
## of the long-run correlation matrix with exact DWT of Fay et al (2009). 
## 
## 	INPUT	 d	 kx1 long-range dependence parameters
##               x	 Data (nxk vector)
##		 filter	 Wavelet filter
##               psih	 List containing psih$fct, the Fourier transform of the 
##			 wavelet mother at values psih$grid
##               LU	 Bivariate vector (optional) containing 
##			 L, the lowest resolution in wavelet decomposition
##               	 U, the maximal resolution in wavelet decomposition
##
##	OUTPUT		Wavelet Whittle criterion
##
##                                           		Achard & Gannaz (2014)
##_________________________________________________________________________________


if(is.matrix(x)){
	N <- dim(x)[1]
	k <- dim(x)[2]
}else{
	N <- length(x)
	k <- 1
}
x <- as.matrix(x,N,k)

## Wavelet decomposition
xwav <- matrix(0,N,k)
for(j in 1:k){
	xx <- x[,j]		   
	resw<-DWTexact(xx,filter)
	xwav_temp <- resw$dwt
	index<-resw$indmaxband
	Jmax <- resw$Jmax
	xwav[1:index[Jmax],j] <- xwav_temp
}

## we free some memory
new_xwav <- matrix(0,min(index[Jmax],N),k)
if(index[Jmax]<N){
	new_xwav[(1:(index[Jmax])),] <- xwav[(1:(index[Jmax])),]
}
xwav <- new_xwav
index <- c(0,index)

## Computation of psih for the computation of K in the paper
res_psi <- psi_hat_exact(filter,Jmax)
psih <- res_psi$psih
grid_K <- res_psi$grid


return(mww_wav_cov_eval(d,xwav,index,psih,grid_K,LU))
}
