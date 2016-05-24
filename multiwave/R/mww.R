mww <- function(x,filter,LU=NULL){

## Computes the multivariate wavelet Whittle estimation for the of 
## the long-memory parameter d and the long-run covariance matrix. 
## 
## 	INPUT	x	Data (n*k vector)
##		filter	Wavelet filter
##              LU	bivariate vector (optional) containing 
##			L, the lowest resolution in wavelet decomposition
##               	U, the maximal resolution in wavelet decomposition
##
##	OUTPUT  d	Long-range parameter estimation
##		cov	Long-run covariance matrix estimation
##
##						Achard & Gannaz (2014)
##______________________________________________________________________________


k <- 1
my_method <- 'Brent'
lower <- -10
upper <- 10
if(is.matrix(x)){ 
	k <- dim(x)[2] 
	my_method <- 'Nelder-Mead'
	lower <- -Inf
	upper <- Inf
}

md <- optim(rep(0,k),mww_eval,x=x,filter=filter,LU=LU,method=my_method,lower=lower,upper=upper)$par
mg <- mww_cov_eval(md,x,filter,LU)

list(d=md,cov=mg)

}
