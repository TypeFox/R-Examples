mww_wav <- function(xwav,index,psih,grid_K,LU=NULL){

## Computes the multivariate Fourier Whittle estimation for the of 
## the long-memory parameter d and the long-run covariance matrix. 
## 
## 	INPUT	 xwav	Vector of wavelet coefficients
##               index	Vector containing the indexes of xwav 
##			where the coefficients change of scale
##               psih	List containing psih$fct, the Fourier transform of the 
##			wavelet mother at values psih$grid
##               LU	bivariate vector (optional) containing 
##			L, the lowest resolution in wavelet decomposition
##               	U, the maximal resolution in wavelet decomposition
##
##	OUTPUT   d	Long-range parameter estimation
##		 cov	Long-run covariance matrix estimation
##
##                                           		Achard & Gannaz (2014)
##_________________________________________________________________________________


k <- 1
my_method<-'Brent'
lower<- -10
upper<-10
if(is.matrix(xwav)){ 
k <- dim(xwav)[2] 
my_method='Nelder-Mead'
lower<- -Inf
upper<- Inf
}

md <- optim(rep(0,k),mww_wav_eval,xwav=xwav,index=index,LU=LU,method=my_method,lower=lower,upper=upper)$par
mg <- mww_wav_cov_eval(md,xwav,index,psih,grid_K,LU)

list(d=md,cov=mg)

}
