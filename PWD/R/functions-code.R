##################################################
### Power Weighted Densities (PWD) Functions

library(Rcpp)
library(RcppArmadillo)



##' Estimates PWD Parameter alpha by Maximum Marginal Predictive Likelihood
##'
##' This is the main function of the package.  It takes as inputs the 
##' time series data as response, as well as a predictor matrix, excluding 
##' the intercept column, and other settings. It returns as outputs a 
##' scalar representing the value of alpha which maximizes the marginal 
##' predictive likelihood of the data given the grid of alpha values 
##' considered.
##' @name alphahat_LR_one_Rcpp
##' @param y T-length time series vector.  y[1] represents the beginning of
##' the time eries.
##' @param X [T x p] dimensional matrix of covariates. This should not 
##' include the intercept column. If X is FALSE, intercept model is run.
##' @param alpha.grid Grid of alpha values over which to compute the marginal
##' predictive likelihood.
##' @param init integer representing the time point to begin computing marginal
##' predictive likelihoods.
##' @param plotting If TRUE, plot the marginal predictive distribution of alpha.
##' @return Return a scalar value representing the value of alpha which maximizes
##' the marginal predictive likelihood of the data over the grid of alpha values
##' considered.
##' @examples 
##' set.seed(12)
##' N=80
##' err = rnorm(N)
##' X = 1:N
##' slopes = c(rep(1.5,40),rep(2,N-40))
##' y = rep(5,N) + slopes*X + err
##' init=6                
##' alpha.grid = seq(.75,1,length.out=40)
##' alphahat=alphahat_LR_one_Rcpp(y=y,X=X,alpha.grid=alpha.grid,init=init,plotting=TRUE)
##' alpha1 = 1.0
##' coeffs1 = bhat.func(y,X,alpha1)
##' alpha2 = alphahat
##' coeffs2 = bhat.func(y,X,alpha2)
##' plot(x=X,y=y)
##' abline(a=coeffs2[1],b=coeffs2[2],lty=2,col="red")
##' abline(a=coeffs1[1],b=coeffs1[2],lty=2)
##' legend("right", legend=c("OLS","PWD"), col=c(1,2), lty=c(2,2), lwd=c(1,1))
alphahat_LR_one_Rcpp = function(y, X=FALSE, alpha.grid=seq(.65,1,length.out=150), init=2, plotting=TRUE){
	# Inputs: 
	# y: time series vector [y1 ... yT]
	# X: [T x p] matrix or data.frame of covariates, not including the intercept column which will be added.
	#	    If X=FALSE, no covariates - intercept model only. The default is that no covariates are included.
	# alpha.grid: Grid of alpha values over which to compute the marginal predictive likelihood.
	# init: time point to begin computing marginal predictive likelihoods.
	# plotting: if TRUE, plot marginal predictive distribution of alpha. 

	# Output: 
	# alpha_hat: scalar point estimator for alpha hat via maximum marginal predictive likelihood
	#		over range of alpha values considered.
	
	sum.a.min = sum(alpha.grid[1]^(0:(init-1)))
	if(!(is.matrix(X)|is.data.frame(X)|is.numeric(X))){p=1}
	if(is.numeric(X)){p=2}
	if((is.matrix(X))|(is.data.frame(X))){p=ncol(X)+1}
	X=as.matrix(X)
	if((sum.a.min-p-1)<=0){
		print("You need to increase the bottom end of alpha.grid, or init, or both - negative DF's.")
		}
	grid.range = alpha.grid
	nsim = length(grid.range)
	full = rep(NA,nsim)
	for(i in 1:nsim){full[i] = loglik.norm.LR.Rcpp(y,X=X,alpha=grid.range[i],init=init)}        
	full[(abs(full)==Inf)|is.na(full)] = min(full) - 1000000000   # defaulting NA/Inf to trivially small log likelihoods
	full = exp(full-max(full))
	full = full/sum(full)
	if(plotting==TRUE){
		ttl = "Marginal Distribution of Alpha"
		plot(x=grid.range,y=full,type="l",col="blue",main=ttl, ylab="Density")
		} 	
    alpha_hat = grid.range[full==max(full)]
	if(alpha_hat==min(grid.range)){
		print("Warning - parameter is minimum of grid - might want to expand grid (and increase init?)")
		}
    alpha_hat
	}



	
	

##' Compute Marginal Predictive Loglikelihood of Data Given alpha
##'
##' This function computes the marginal predictive loglikelihood of the 
##' observed data given a particular value of alpha, and the time
##' point to begin computing marginal predictive loglikelihoods.
##' @name loglik.norm.LR.Rcpp
##' @param y T-length time series vector.  y[1] represents the beginning of
##' the time eries.
##' @param X [T x p] dimensional matrix of covariates. This should not 
##' include the intercept column. If X is FALSE, intercept model is run.
##' @param alpha PWD parameter we are calculating the marginal predictive
##' loglikelihood for.
##' @param init integer representing the time point to begin computing marginal
##' predictive loglikelihood.
##' @return Return a scalar value representing the marginal predictive 
##' loglikelihood of the data given alpha.
##' @examples 
##' N=80
##' err = rnorm(N)
##' X = 1:N
##' slopes = c(rep(1.5,40),rep(2,N-40))
##' y = rep(5,N) + slopes*X + err                  
##' init=4
##' alpha.grid = seq(.65,1,length.out=40)
##' i=40
##' loglik.norm.LR.Rcpp(y,X=X,alpha=alpha.grid[i],init=init)
loglik.norm.LR.Rcpp = function(y,X=FALSE,alpha,init){

	# Inputs: 
	# y: time series vector [y1 ... yT]
	# 		y[1] should be y_1 (i.e., the time point furthest into the past).
	# X: matrix or data.frame of covariates, not including the intercept column which will be added.
	#	   	If X=FALSE, no covariates - intercept model only.
	# alpha: PWD parameter we are calculating the marginal predictive loglikelihood for.
	# init: time point to begin computing marginal predictive loglikelihoods.
	
	# Output: 
	# LL: scalar-valued marginal predictive loglikelihood, P(y|X,alpha)
	
	n = length(y)
	output = 0
	icept=rep(1,length(y))
	
	# Add intercept vector:
	if(!(is.matrix(X)|is.data.frame(X)|is.numeric(X))){X_aug = as.matrix(icept,ncol=1)}
	if(is.matrix(X)|is.data.frame(X)|is.numeric(X)){X_aug = cbind(icept,X)}
		
	LL = logliknormLR(yy=y,XX_aug=X_aug,alpha=alpha,init=init)	
	LL
	}
	

##' Compute PWD Regression Coefficients Given alpha
##'
##' This function computes PWD regression coefficients for response
##' y and predictors X given a particular value of alpha.
##' @param y T-length time series vector.  y[1] represents the beginning of
##' the time eries.
##' @param X [T x p] dimensional matrix of covariates. This should not 
##' include the intercept column. If X is FALSE, intercept model is run.
##' @param alpha PWD parameter we are calculating the marginal predictive
##' loglikelihood for.
##' @return (p+1)-length vector representing the regression coefficients
##' associated with a PWD regression of y upon X given PWD parameter alpha.
##' @examples 
##' set.seed(12)
##' N=80
##' err = rnorm(N)
##' X = 1:N
##' slopes = c(rep(1.5,40),rep(2,N-40))
##' y = rep(5,N) + slopes*X + err                  
##' alpha1 = 1.0
##' coeffs1 = bhat.func(y,X,alpha1)
##' alpha2 = .9
##' coeffs2 = bhat.func(y,X,alpha2)
##' plot(x=X,y=y)
##' abline(a=coeffs2[1],b=coeffs2[2],lty=2,col="red")
##' abline(a=coeffs1[1],b=coeffs1[2],lty=2)
##' legend("right", legend=c("OLS","PWD"), col=c(1,2), lty=c(2,2), lwd=c(1,1))
bhat.func = function(y, X, alpha){
	# Inputs: 
	# y: single vector time series [y1 ... yT]
	# X: matrix or data.frame of covariates, not including the intercept column which will be added.
	#	   	If X=FALSE, no covariates - intercept model only.
	# alpha: Estimated alpha value for time series.
	
	# Output: 
	# bhat: vector of PWD regression coefficient estimates.
	
	T = length(y)
	icept=rep(1,T)
	if(!(is.matrix(X)|is.data.frame(X)|is.numeric(X))){X = rep(1,T)}
	bhat = as.numeric(summary(lm(y~X, weights=alpha^((T-1):0)))$coefficients[,"Estimate"])   # equivalent to WLS.
	bhat
	}




##' Fast Computation of Marginal Predictive Loglikelihood
##'
##' Helper function which takes as input a vector-valued response, y, 
##' a predictor matrix, X, a particular value of alpha, and init. It
##' returns as an output the marginal predictive loglikelihood of the 
##' data given that value of alpha.
##' @name logliknormLR
##' @param yy T-length time series vector.  y[1] represents the beginning of
##' the time eries.
##' @param XX_aug [T x (p+1)] dimensional matrix of covariates. This will
##' include the intercept column. 
##' @param alpha PWD parameter we are calculating the marginal predictive
##' loglikelihood for.
##' @param init integer representing the time point to begin computing marginal
##' predictive loglikelihood.
##' @return Return a scalar value representing the marginal predictive 
##' loglikelihood of the data given alpha.
logliknormLR <- function(yy, XX_aug, alpha, init) {
    .Call('PWD_logliknormLR', PACKAGE = 'PWD', yy, XX_aug, alpha, init)
}

	
	
	
	
	

