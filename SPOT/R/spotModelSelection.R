###################################################################################################
#' Model selection and error estimation
#' 
#' This is a number of functions for model selection and error estimation: \cr
#' rsq: R-Squared, can be adjusted with a complexity measure \cr
#' sse: Sum of Squared Errors \cr
#' sae: Sum of Absolute Errors \cr
#' mse: Mean Squared Errors \cr
#' rmse: Root Mean Squared Errors \cr
#' aic: Akaike Information Criterion, with bias adjustment for small sample sizes \cr
#' scaled: this refers to an approach that scales the data before error calculation, see Keijzer (2004).
#'
#' @name spotSelectionCriteria
#' @usage spotSelectionRsq(yi,fi)
#' spotSelectionAdjustedRsq(yi,fi,p)
#' spotSelectionAic(yi,fi,p)
#' spotSelectionSse(yi,fi)
#' spotSelectionSae(yi,fi)
#' spotSelectionMse(yi,fi)
#' spotSelectionMae(yi,fi)
#' spotSelectionRmse(yi,fi)
#' spotSelectionScaledSse(yi,fi)
#' spotSelectionScaledMse(yi,fi)
#' spotSelectionScaledRmse(yi,fi); 
#'
#' @param yi sampled values vector
#' @param fi predicted values vector
#' @param p complexity measure or number of regressors
#'
#' @return a scalar value is returned
#'
#' @export spotSelectionRsq spotSelectionAdjustedRsq spotSelectionAic spotSelectionSse spotSelectionSae spotSelectionMse spotSelectionMae spotSelectionRmse spotSelectionScaledSse spotSelectionScaledMse spotSelectionScaledRmse
#' @aliases spotSelectionRsq spotSelectionAdjustedRsq spotSelectionAic spotSelectionSse spotSelectionSae spotSelectionMse spotSelectionMae spotSelectionRmse spotSelectionScaledSse spotSelectionScaledMse spotSelectionScaledRmse
#' @references - Maarten Keijzer. 2004. Scaled Symbolic Regression. Genetic Programming and Evolvable Machines 5 (3): 259-269. \cr
#' -  Hirotugu Akaike. 1974. A new look at the statistical model identification. IEEE Transactions on Automatic Control 19 (6): 716-723.
###################################################################################################
spotSelectionRsq <- function(yi,fi){
	ybar<-mean(yi)
	SStot<-sum((yi-ybar)^2)
	SSerr<-sum((yi-fi)^2)
	-(1-(SSerr/SStot))
}

# adjusted r-squared 
spotSelectionAdjustedRsq <- function(yi,fi,p){ #p is number of regressors
	nn <- length(yi) #nn is sample size
	if(nn <= (p+1)){return(Inf) #does not work properly if number p is too damn high, TODO
	}else{return(-(-spotSelectionRsq(yi,fi)-(p/(nn-1)))*(nn-1)/(nn-p-1))}
}

# sum of squared errors
spotSelectionSse <- function(yi,fi){ 	sum((yi-fi)^2) }

# sum of absolute errors
spotSelectionSae <- function(yi,fi){ 	sum(abs(yi-fi)) }

# mean squared error
spotSelectionMse <- function(yi,fi){ 	(1/length(yi)) * spotSelectionSse(yi,fi) }

# mean absolute error
spotSelectionMae <- function(yi,fi){ 	(1/length(yi)) * spotSelectionSae(yi,fi) }

# root mean squared error
spotSelectionRmse <- function(yi,fi){	sqrt(spotSelectionMse(yi,fi)) }


# scaled sum of squared errors
spotSelectionScaledSse <- function(yi,fi){
	if(length(fi)==1)fi<-rep(fi,length(yi))
	varfi<-var(fi)
	if(varfi==0){
		return(sqrt(mean((yi-mean(yi))^2))) #fi is a constant, constant is shifted to mean of yi
	}else{
		b=cov(yi,fi)/var(fi)
		a=mean(yi)-b*mean(fi)
		return(sum((yi-(a+b*fi))^2))
	}
}

#http://www2.cs.uidaho.edu/~cs472_572/f11/scaledsymbolicRegression.pdf
# scaled mean squared error
spotSelectionScaledMse <- function(yi,fi){	(1/length(yi)) * spotSelectionScaledSse(yi,fi) }
spotSelectionScaledRmse <- function(yi,fi){	sqrt(spotSelectionScaledMse(yi,fi)) }


#this is aic with bias adjustment for small sample sizes

# akaike information criterion
spotSelectionAic <- function(yi,fi,p){ #p is number of regressors
	nn <- length(yi) #nn is sample size
	SSerr<-spotSelectionSse(yi,fi)
	#return(nn*log(SSerr/nn)+2*p) #no bias adjustment for small sample size
	if(nn <= (p+1)){return(Inf) #does not work properly if number p is too damn high, TODO
	}else{return(nn*log(SSerr/nn)+2*p+ (2*p*(p+1))/(nn-p-1))} #smaller is better
}