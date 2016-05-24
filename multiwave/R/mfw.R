mfw <- function(x,m){

## Computes the multivariate Fourier Whittle estimation for the of 
## the long-memory parameter d and the long-run covariance matrix. 
## 
## 	INPUT	x	Data (n*k matrix)
##		m	Truncation number in Fourier frequencies
##			
##	OUTPUT  d	Long-range parameter estimation
##		cov	Long-run covariance matrix estimation
##
##					based on the paper of Katsumi Shimotsu, 2007
##					Achard & Gannaz (2014)
##__________________________________________________________________________________

k <- 1
my_method<-'Brent'
lower<- -10
upper<-10
if(is.matrix(x)){ 
k <- dim(x)[2] 
my_method='Nelder-Mead'
lower<- -Inf
upper<- Inf
}

md <- optim(rep(0,k),mfw_eval,x=x,m=m,method=my_method,lower=lower,upper=upper)$par
mg <- mfw_cov_eval(md,x,m)

list(d=md,cov=mg)

}
