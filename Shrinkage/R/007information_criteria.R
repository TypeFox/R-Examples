
# Corey M. Yanofsky, February 18, 2009
# M. Padilla, Oct, 2013 (small modifications)
########################################################################################################### 

# AICc-, BIC-, and Bayes-factor-based statistics

logBF2posprob <- function(logBayesFactor, prior.odds = 1){
# This function takes evidence on the Jeffreys scale (the argument logBayesFactor) and transforms it onto the probability scale.
# If logBayesFactor comes from BIC.one.covariate(...)*.5 or AICc.one.covariate(...)*.5, then the result is the approximate 
# probability (or Akaike weight, for AICc) in favor of the smaller model.

	1 - (1 / (1 + prior.odds*exp(logBayesFactor)))
}
setGeneric("BIC.one.covariate", function(x,y,...) standardGeneric("BIC.one.covariate"))
setGeneric("BF.one.covariate", function(x,y,...) standardGeneric("BF.one.covariate"))
setGeneric("AICc.one.covariate", function(x,y,...) standardGeneric("AICc.one.covariate"))
#---------------
setMethod("BIC.one.covariate", signature(x = "numeric", y = "missing"), function(x, y, design_matrix, covariate_vector, include_intercept = TRUE) 
{
# The returned value is the difference of the BICs for a base model versus one with one additional covariate.
# Positive returned values favor the smaller model.
# Use logBF2posprob(.5*BIC.one.covariate(...)) to get the approximate Bayes factor for the base model. The result is already normalized, so 
# 1 - logBF2posprob(.5*BIC.one.covariate(...)) is the approximate Bayes factor for the expanded model.
#
# x - the data vector for the dependent variable 
#
# design matrix - the design matrix for the smaller model; can be a data frame. Note that the name "x" is reserved for the dependent variable.
#
# covariate_vector - the vector of covariate values for the parameter in the expanded model
#
# include_intercept - a logical variable for including an intercept in the model; defaults to TRUE 
#
# if neither design_matrix nor covariate_vector is given, then the BIC for zero vs. non-zero mean is calculated (include_intercept is ignored)
# if design_matrix is missing, covariate_vector is given, and include_intercept == TRUE,
#	then the BIC for mu = mu_0 vs. mu = mu_0 + covariate_vector*beta is calculated
# if design_matrix is missing, covariate_vector is given, and include_intercept == FALSE,
#	then the BIC for mu = 0 vs. mu = covariate_vector*beta is calculated
# it is an error to give design_matrix and not covariate_vector



	stopifnot(!(missing(covariate_vector) && !missing(design_matrix)))
	n_data <- sum(!is.na(x))
	if(n_data == 0)
		return(0)
	if(!missing(covariate_vector))
	{
		if(missing(design_matrix))
		{
			SSR_0 <- (n_data-1)*var(na.omit(x))
			if(include_intercept)
				data_1 <- data.frame(cbind(x,rep(1,length(x)),covariate_vector))
			else
				data_1 <- data.frame(cbind(x,covariate_vector))
		}
		else
		{	
			data_0 <- data.frame(cbind(x,design_matrix))
			if(include_intercept)
			{
				data_1 <- data.frame(cbind(x,rep(1,length(x)),design_matrix,covariate_vector))
				fit_0 <- lm(x ~ ., data_0, na.action = na.omit)
			}
			else
			{
				data_1 <- data.frame(cbind(x,design_matrix,covariate_vector))
				fit_0 <- lm(x ~ . + 0, data_0, na.action = na.omit)
			}
			SSR_0 <- sum(residuals(fit_0)^2)
			
		}
		fit_1 <- lm(x ~ . + 0, data_1, na.action = na.omit)
		SSR_1 <- sum(residuals(fit_1)^2)
	}
	else if(missing(design_matrix))
	{	# then test for non-zero mean
		x <- na.omit(x)
		SSR_0 <- sum(x^2)
		SSR_1 <- (n_data-1)*var(x)
	}
	
	return_val = (n_data*(log(SSR_1) - log(SSR_0)) + log(n_data))
	
	if(is.na(return_val))
	{
			warning("not enough data to do model selection")
			return(0)
	}
	
	return(return_val)
		
})

setMethod("BIC.one.covariate", signature(x = "numeric", y = "numeric"), function(x, y, x_design_matrix, y_design_matrix) 
{
# returns the BIC for testing if x and y have the same means, conditional on the covariates
#
# x - the data vector for group 1
#
# y - the data vector for group 2
#
# x_design matrix - the design matrix for group 1
#
# y_design matrix - the design matrix for group 2
#
# if design matrices are given, they must have the same number of covariates (columns)
# it is an error to give only one design matrix

	stopifnot(missing(x_design_matrix) && missing(y_design_matrix) || (!missing(x_design_matrix) && !missing(y_design_matrix)))
	covariate_vector <- c(rep(0,length(x)),rep(1,length(y)))
	x <- c(x,y)
	if(!missing(x_design_matrix))
	{
		design_matrix = rbind(x_design_matrix,y_design_matrix)
		BIC.one.covariate(x = x, design_matrix = design_matrix, covariate_vector = covariate_vector)
	}
	else
		BIC.one.covariate(x = x, covariate_vector = covariate_vector)

})
#---------------
setMethod("BF.one.covariate", signature(x = "numeric", y = "missing"), function(x, y, prior_mean = 0, design_matrix, covariate_vector, include_intercept = TRUE) 
{
# Returns the logarithm Bayes factor (i.e., evidence on the Jeffreys scale) for one additional parameter. Positive values favor the smaller model.
# The prior in the expanded model is the unit-information prior.
# Use logBF2posprob(BF.one.covariate(...)) to get the Bayes factor for the base model. The result is already normalized, so 
# 1 - logBF2posprob(BF.one.covariate(...)) is the Bayes factor for the expanded model.
#
# x - the data vector for the dependent variable 
#
# prior_mean - the prior mean for the unit-information prior on the new parameter; default value is 0
#
# design matrix - the design matrix for the smaller model; can be a data frame. Note that the name "x" is reserved for the dependent variable.
#
# covariate_vector - the vector of covariate values for the parameter in the expanded model
#
# include_intercept - a logical variable for including an intercept in the model; defaults to TRUE 
#
# if neither design_matrix nor covariate_vector is given, then the BF for zero vs. non-zero mean is calculated (include_intercept is ignored)
# if design_matrix is missing, covariate_vector is given, and include_intercept == TRUE,
#	then the BF for mu = mu_0 vs. mu = mu_0 + covariate_vector*beta is calculated
# if design_matrix is missing, covariate_vector is given, and include_intercept == FALSE,
#	then the BF for mu = 0 vs. mu = covariate_vector*beta is calculated
# it is an error to give design_matrix and not covariate_vector

	stopifnot(!(missing(covariate_vector) && !missing(design_matrix)))
	n_data <- sum(!is.na(x))
	if(n_data == 0)
		return(0)
	if(!missing(covariate_vector))
	{
		if(missing(design_matrix))
		{
			n_param <- as.numeric(include_intercept)
			SSR_0 <- (n_data-1)*var(na.omit(x))
			if(include_intercept)
				data_1 <- data.frame(rbind(cbind(x,rep(1,length(x)),covariate_vector), c(prior_mean, 0, 1)))
			else
				data_1 <- data.frame(rbind(cbind(x,covariate_vector), c(prior_mean, 1)))
			complexity_0 <- .5*log(n_data)
		}
		else
		{	
			n_param <- ncol(design_matrix) + as.numeric(include_intercept)
			data_0 <- data.frame(cbind(x,design_matrix))
			if(include_intercept)
			{
				data_1 <- data.frame(rbind(cbind(x,rep(1,length(x)),design_matrix,covariate_vector), c(prior_mean, rep(0, n_param), 1)))
				fit_0 <- lm(x ~ ., data_0, na.action = na.omit, qr = TRUE)
			}
			else
			{
				data_1 <- data.frame(rbind(cbind(x,design_matrix,covariate_vector), c(prior_mean, rep(0, n_param), 1)))
				fit_0 <- lm(x ~ . + 0, data_0, na.action = na.omit, qr = TRUE)
			}
			SSR_0 <- sum(residuals(fit_0)^2)
			complexity_0 = sum(log(abs(diag(qr.R(fit_0$qr)))))
		}
		fit_1 <- lm(x ~ . + 0, data_1, na.action = na.omit)
		SSR_1 <- sum(residuals(fit_1)^2)
		# sum(log(abs(diag(qr.R(fit$qr))))) = .5.*log(det(t(X)%*%X)), where X is the matrix of covariates
		complexity_penalty <- sum(log(abs(diag(qr.R(fit_1$qr))))) - complexity_0
	}
	else if(missing(design_matrix))
	{	# then test for non-zero mean
		n_param <- 0
		x <- na.omit(x)
		SSR_0 <- sum(x^2)
		SSR_1 <- n_data*var(c(prior_mean, x))
		complexity_penalty <- .5*log(n_data+1)
	}
	
	complexity_penalty + .5*(n_data - n_param)*(log(SSR_1) - log(SSR_0))
})
setMethod("BF.one.covariate", signature(x = "numeric", y = "numeric"), function(x, y, prior_mean = 0, x_design_matrix, y_design_matrix) 
{
# Returns the logarithm of the Bayes factor for testing if x and y have the same means, conditional on the covariates.
#
# x - the data vector for group 1
#
# y - the data vector for group 2
#
# prior_mean - the prior mean for the unit-information prior on the difference of the means; default value is 0
#
# x_design matrix - the design matrix for group 1
#
# y_design matrix - the design matrix for group 2
#
# if design matrices are given, they must have the same number of covariates (columns)
# it is an error to give only one design matrix

	stopifnot(missing(x_design_matrix) && missing(y_design_matrix) || (!missing(x_design_matrix) && !missing(y_design_matrix)))
	if(all(is.na(x)) || all(is.na(y)))
		return(0)
	covariate_vector <- c(rep(0,length(x)),rep(1,length(y)))
	x <- c(x,y)
	if(!missing(x_design_matrix))
	{
		design_matrix = rbind(x_design_matrix,y_design_matrix)
		BF.one.covariate(x = x, prior_mean = prior_mean, design_matrix = design_matrix, covariate_vector = covariate_vector)
	}
	else
		BF.one.covariate(x = x, prior_mean = prior_mean, covariate_vector = covariate_vector)

})

#---------------
setMethod("AICc.one.covariate", signature(x = "numeric", y = "missing"), function(x, y, design_matrix, covariate_vector, include_intercept = TRUE) 
{
# The returned value is the difference of the AICc values for a base model versus one with one additional covariate.
# Positive returned values favor the smaller model.
# Use logBF2posprob(.5*AICc.one.covariate(...)) to get the Akaike weight for the base model. The result is already normalized, so 
# 1 - logBF2posprob(.5*AICc.one.covariate(...)) is the Akaike weight for the expanded model.
#
# x - the data vector for the dependent variable 
#
# design matrix - the design matrix for the smaller model; can be a data frame. Note that the name "x" is reserved for the dependent variable.
#
# covariate_vector - the vector of covariate values for the parameter in the expanded model
#
# include_intercept - a logical variable for including an intercept in the model; defaults to TRUE 
#
# if neither design_matrix nor covariate_vector is given, then the AICc for zero vs. non-zero mean is calculated (include_intercept is ignored)
# if design_matrix is missing, covariate_vector is given, and include_intercept == TRUE,
#	then the BF for mu = mu_0 vs. mu = mu_0 + covariate_vector*beta is calculated
# if design_matrix is missing, covariate_vector is given, and include_intercept == FALSE,
#	then the BF for mu = 0 vs. mu = covariate_vector*beta is calculated
# it is an error to give design_matrix and not covariate_vector

	stopifnot(!(missing(covariate_vector) && !missing(design_matrix)))
	n_data <- sum(!is.na(x))
	if(n_data == 0)
		return(0)
	if(!missing(covariate_vector))
	{
		if(missing(design_matrix))
		{
			SSR_0 <- (n_data-1)*var(na.omit(x))
			K_0 <- 2 # parameters - mean and sd
			if(include_intercept)
				data_1 <- data.frame(cbind(x,rep(1,length(x)),covariate_vector))
			else
				data_1 <- data.frame(cbind(x,covariate_vector))
		}
		else
		{	
			data_0 <- data.frame(cbind(x,design_matrix))
			if(include_intercept)
			{
				data_1 <- data.frame(cbind(x,rep(1,length(x)),design_matrix,covariate_vector))
				fit_0 <- lm(x ~ ., data_0, na.action = na.omit)
				K_0 <- ncol(data_0) + 1 # parameters - one per covariate column of data_0 and intercept and sd
			}
			else
			{
				data_1 <- data.frame(cbind(x,design_matrix,covariate_vector))
				fit_0 <- lm(x ~ . + 0, data_0, na.action = na.omit)
				K_0 <- ncol(data_0) # parameters - one per covariate column of data_0 and sd
			}
			SSR_0 <- sum(residuals(fit_0)^2)
			
		}
		fit_1 <- lm(x ~ . + 0, data_1, na.action = na.omit)
		K_1 <- ncol(data_1) # parameters - one per covariate column of data_1 and sd
		SSR_1 <- sum(residuals(fit_1)^2)
	}
	else if(missing(design_matrix))
	{	# then test for non-zero mean
		x <- na.omit(x)
		SSR_0 <- sum(x^2)
		SSR_1 <- (n_data-1)*var(x)
		K_0 <- 1 # parameters - sd
		K_1 <- 2 # parameters - mean and sd
	}
	if((n_data - K_1 - 1) <= 0) 
	{
		warning("not enough data to do model selection")
		return(0)
	}
	
	n_data*(log(SSR_1) - log(SSR_0)) + 2*K_1 + (2*K_1*(K_1 + 1)/(n_data - K_1 - 1)) - 2*K_0 - (2*K_0*(K_0 + 1)/(n_data - K_0 - 1))
})
setMethod("AICc.one.covariate", signature(x = "numeric", y = "numeric"), function(x, y, x_design_matrix, y_design_matrix) 
{
# returns the AICc for testing if x and y have the same means, conditional on the covariates
#
# x - the data vector for group 1
#
# y - the data vector for group 2
#
# x_design matrix - the design matrix for group 1
#
# y_design matrix - the design matrix for group 2
#
# if design matrices are given, they must have the same number of covariates (columns)
# it is an error to give only one design matrix

	stopifnot(missing(x_design_matrix) && missing(y_design_matrix) || (!missing(x_design_matrix) && !missing(y_design_matrix)))
	covariate_vector <- c(rep(0,length(x)),rep(1,length(y)))
	x <- c(x,y)
	if(!missing(x_design_matrix))
	{
		design_matrix = rbind(x_design_matrix,y_design_matrix)
		AICc.one.covariate(x = x, design_matrix = design_matrix, covariate_vector = covariate_vector)
	}
	else
		AICc.one.covariate(x = x, covariate_vector = covariate_vector)

})

#---------------
BICfdr <- function(x, y=NULL,prior.odds = 1, ...) {
	if(is(x,"matrix")){x<-as.numeric(x);x<-MakeNames(x)}
	if(is(y,"matrix")){y<-as.numeric(y);y<-MakeNames(y)}
	if(length(y)==0){
		zo<-nNumeric(logBF2posprob(.5*BIC.one.covariate(x = x,...),prior.odds = 1))	
	}
	else{
		zo<-nNumeric(logBF2posprob(.5*BIC.one.covariate(x = x,y=y,...),prior.odds = 1))	
	}
}
BFfdr <- function(x, y=NULL,prior.odds = 1, ...) {
	if(is(x,"matrix")){x<-as.numeric(x);x<-MakeNames(x)}
	if(is(y,"matrix")){y<-as.numeric(y);y<-MakeNames(y)}
	if(length(y)==0){
		zo<-nNumeric(logBF2posprob(BF.one.covariate(x = x,...),prior.odds = 1))	
	}
	else{
		zo<-nNumeric(logBF2posprob(BF.one.covariate(x = x,y=y,...),prior.odds = 1))	
	}
}
AICcfdr <- function(x, y=NULL,prior.odds = 1, ...) {
	if(is(x,"matrix")){x<-as.numeric(x);x<-MakeNames(x)}
	if(is(y,"matrix")){y<-as.numeric(y);y<-MakeNames(y)}
	if(length(y)==0){
		zo<-nNumeric(logBF2posprob(.5*AICc.one.covariate(x = x,...),prior.odds = 1))	
	}
	else{
		zo<-nNumeric(logBF2posprob(.5*AICc.one.covariate(x = x,y=y,...),prior.odds = 1))	
	}
}
#----------------
