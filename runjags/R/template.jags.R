#' @title Generate a generalised linear mixed model (GLMM) specification in JAGS
#' @name template.jags
#' @aliases template.jags template.JAGS
#' @export

#' @description
#' Use an lme4 style syntax to create a JAGS model representation of a GLMM, including all data, initial values and monitor specifications required to run the model using \code{\link{run.jags}}.

#' @details
#' This function is designed to allow new users to MCMC to create relatively simple GLMM models in JAGS using an lme4-style formula interface.  Examining the template created by this function is a good way to learn about how the BUGS language is structured, as well as the options provided by the runjags package.  After generating the template model, the user is encouraged to examine the model file and make whatever changes are necessary before running the model using `run.jags'.  You can also run the models with no changes and compapre the results to those obtained through more standard model fitting approaches to learn more about how the differently presented sets of inference relate to each other.  Note that the effect of the reference level for factors is explicitly given as 0 in output from runjags.  For more about the BUGS language, see Lunn et al (2012).

#' @keywords models

#' @return
#' The filename of the created model template.

#' @seealso
#' \code{\link{run.jags}} to run the model, \code{\link{add.summary}} for details of summary statistics available from the fitted model, and \code{\link{runjags-class}} for details of how to extract information such as residuals and the fitted values.

#' @references 
#' Lunn D, Jackson C, Best N, Thomas A, Spiegelhalter D (2012). The BUGS book: A practical introduction to Bayesian analysis. CRC press.

#' @examples
#' \dontrun{
#' # Create a simple linear model and compare the results to LM:
#' 
#' # This is based on the example in ?lm:
#' ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
#' trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
#' group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
#' weight <- c(ctl, trt)
#' D9 <- data.frame(weight, group)
#' lm.D9 <- lm(weight ~ group, data=D9)
#' 
#' # The JAGS equivalent:
#' model <- template.jags(weight ~ group, D9, n.chains=2, 
#' family='gaussian')
#' JAGS.D9 <- run.jags(model)
#' summary(JAGS.D9)
#' summary(lm.D9)
#' # Note that lm reports sigma and JAGS the precision - to
#' # make them more comparable we could use a mutate function:
#' JAGS.D9 <- run.jags(model, mutate=list(prec2sd, 'precision'))
#' summary(JAGS.D9)
#' summary(lm.D9)
#' # Compare the estimated residuals:
#' plot(residuals(lm.D9), residuals(JAGS.D9, output='mean'))
#' 
#' # For more examples see:
#' vignette('quickjags', package='runjags')
#' }

#' @param formula a formula representation of the desired model, using lme4 style syntax.  Two-way interactions for all variables are permitted, as are random intercepts.

#' @param data a data frame containing the variables specified in formula.  This must be specified.

#' @param file the filename of the model to output.  This will be over-written if it exists.

#' @param family a character string representing the response distribution - options are:  'gaussian', 'binomial', 'Poisson', 'negative binomial', 'ZIB', 'ZIP', 'ZINB' (the latter denote zero-inflated distributions).

#' @param write.data option to write the data to file with the model.  If the data is very large it may be better not to write this to file, but the same data frame must be given to the subsequent run.jags call that runs the model.

#' @param write.inits option to write the initial values to file with the model.

#' @param precision.prior the prior distribution to be used for precision parameters.

#' @param effect.prior the prior distribution to be used for linear and fixed effect terms, as well as interactions and the intercept.

#' @param n.chains the number of chains to use.

#' @param precision.inits a numeric vector of initial values from which the precision parameters in the model will be randomly chosen.  It is recommended to make these over-dispersed, but if the values are too extreme the model may not compile.

#' @param effect.inits a numeric vector of initial values from which the effect parameters in the model will be randomly chosen.  It is recommended to make these over-dispersed, but if the values are too extreme the model may not compile.

#' @param inits an optional list of named lists to specify initial values for one or more parameters in each chain.  The number of named lists must match n.chains.
NULL

#' @rdname template.jags
template.jags <- function(formula, data, file='JAGSmodel.txt', family='gaussian', write.data=TRUE, write.inits=TRUE, precision.prior='dgamma(0.001, 0.001)', effect.prior='dnorm(0, 10^-6)', n.chains=2, precision.inits=c(0.01,10), effect.inits=c(-1, 1), inits=NULL){
	
	formula <- as.formula(formula)
	data <- as.data.frame(data)
	
	if(length(as.character(formula))!=3)
		stop('Unsupported formula expression')
	
	if(length(precision.inits)<2 || length(effect.inits)<2)
		stop('At least two precision and effect initial values must be provided')
	if(n.chains<2)
		stop('Two or more chains are required for standard methods of assessing convergence')
	
	# Sort out any initial values:
	if(is.null(inits))
		inits <- lapply(1:n.chains, function(x) return(list()))
	if(is.list(inits) && !any(sapply(inits, is.list)))
		inits <- lapply(1:n.chains, function(x) return(inits))
	if(!is.list(inits) || length(inits)!=n.chains || !all(sapply(inits, is.list)))
		stop('If initial values are provided to the template.jags function this must be as a named list of length equal to the number of chains', call.=FALSE)
	passedinits <- inits
	
	# May change this, but only for write.data:
	convertstrings <- FALSE
	
	if(class(family)!='character')
		stop('Invalid family specification - a character string must be supplied')
	possibles <- c('gaussian','binomial','poisson','nb','negative binomial','zib','zip','zinb')
	family <- possibles[pmatch(tolower(family), possibles)]
	if(length(family)!=1 || is.na(family))
		stop('Invalid family specification - consult the help file for the possibilites.  Additional families (and alternative link functions) can be used by manually editing a template file created with a supported family.')
	if(family=='negative binomial')
		family <- 'nb'
	
	zifamily <- FALSE
	if(family=='zib'){
		zifamily <- TRUE
		family <- 'binomial'
	}
	if(family=='zip'){
		zifamily <- TRUE
		family <- 'poisson'
	}
	if(family=='zinb'){
		zifamily <- TRUE
		family <- 'nb'
	}
	
	
  	terms <- gsub('[[:space:]]','', attr(terms(formula), 'term.labels'))
	if(any(grepl('+',terms,fixed=TRUE) | grepl('*',terms,fixed=TRUE)))
		stop('There was a problem parsing the formula - did you forget the parentheses around a random effects term e.g. (1 | random) ?')
	intercept <- attr(terms(formula), 'intercept')
	response <- as.character(formula)[2]

	if(grepl('+',response,fixed=TRUE) || grepl('(',response,fixed=TRUE) || grepl('[',response,fixed=TRUE))
		stop('Unsupported response expression')
	if(any(nchar(gsub('[^:]','',terms))>1))
		stop('Unsupported 3+ way interaction term')
	if(any(grepl('I(',terms,fixed=TRUE)))
		stop('The I( ) construct is not supported - provide variables in a pre-calculated form, or manually edit the template file')
	if(!any(response==names(data)))
		stop('The response variable was not found in the data')
	
	offsets <- gsub('[[:space:]]','',strsplit(as.character(formula)[3],'+',fixed=TRUE)[[1]])
	offsets <- gsub('offset(','', offsets[grepl('offset(', offsets, fixed=TRUE)], fixed=TRUE)
	if(any(grepl('I(',offsets,fixed=TRUE)))
		stop('The I( ) construct is not supported - provide variables in a pre-calculated form, or manually edit the template file')
	if(any(grepl('(',offsets,fixed=TRUE)))
		stop('Functions of offset terms (and parentheses) are not supported - provide variables in a pre-calculated form, or manually edit the template file')
	offsets <- gsub(')','',offsets,fixed=TRUE)
	notfound <- ! offsets %in% names(data)
	if(any(notfound))
		stop(paste('The following offset term(s) was/were not found in the data: ', paste(offsets[notfound], collapse=', ')))	
	
	warnmissing <- FALSE
		
	randoms <- terms[grepl('|',terms,fixed=TRUE)]
	if(any(sapply(strsplit(randoms,'|',fixed=TRUE), length)!=2))
		stop('Random slope terms are not supported by the template function - but you can add these by manually editing the template file created with only random intercept terms')
	if(any(sapply(strsplit(randoms,'|',fixed=TRUE), function(x) return(gsub('[[:space:]]','',x[1])))!='1'))
		stop('Random slope terms are not supported by the template function - but you can add these by manually editing the template file created with only random intercept terms')
	
	madefactors <- character(0)
	if(length(randoms)>0){
		if(grepl('+',randoms,fixed=TRUE) || grepl('/',randoms,fixed=TRUE) || grepl(':',randoms,fixed=TRUE) || grepl('*',randoms,fixed=TRUE) || grepl('-',randoms,fixed=TRUE))
			stop('Unsupported random intercept expression - random effects must take the form (1 | group) - consult the help file for more information')
	
		randoms <- sapply(strsplit(randoms,'|',fixed=TRUE), function(x) return(gsub('[[:space:]]','',x[2])))	
		notfound <- ! randoms %in% names(data)
		if(any(notfound))
			stop(paste('The following random effects term(s) was/were not found in the data: ', paste(randoms[notfound], collapse=', ')))
		for(i in 1:length(randoms)){
			if(class(data[[randoms[i]]])!='factor'){
				data[[randoms[i]]] <- factor(data[[randoms[i]]])
				madefactors <- c(madefactors, randoms[i])
			}
			if(any(is.na(data[[randoms[i]]])))
				warnmissing <- TRUE
		}
	}
	if(length(madefactors)>0)
		warning('One or more random effects terms were coerced into factors')
			
	termvars <- lapply(strsplit(terms[!grepl('|',terms,fixed=TRUE)],':',fixed=TRUE), sort)
	allvars <- unique(unlist(termvars))
	if(any(grepl('I(',allvars,fixed=TRUE)))
		stop('The I( ) construct is not supported - provide variables in a pre-calculated form')
	if(any(grepl('(',allvars,fixed=TRUE)))
		stop('Functions of terms (and parentheses) are not supported - provide variables in a pre-calculated form')
	
	notfound <- ! allvars %in% names(data)
	if(any(notfound))
		stop(paste('The following term(s) is/are not in the data: ', paste(allvars[notfound], collapse=', ')))
	warnmissing <- warnmissing || any(sapply(data[allvars], function(x) return(any(is.na(x)))))
	
	classes <- sapply(data[allvars], class)
	# mode() rather than class() makes integer and double both numeric, but also makes factor numeric
	classes[classes=='integer'] <- 'numeric'
	for(i in which(classes=='character')){
		data[allvars[i]] <- factor(data[allvars[i]])
		madefactors <- c(madefactors, allvars[i])
		classes[i] <- 'factor'
		
		# Not currently allowed:
		stop('Some of the linear and/or fixed effects in the data frame are character variables - please convert these to factors manually, choosing the most appropriate category as the reference level')
	}
	
		
	supported <- classes %in% c('numeric','factor')
	if(any(!supported))
		stop(paste('The following unsupported term classes were found in the data: ', paste(classes[!supported], collapse=', ')))
	
	names(classes) <- allvars
	
	# Check centering of numerics:
	variances <- numeric(0)
	centerwarn <- FALSE
	for(i in which(classes=='numeric')){
		mean <- mean(data[[allvars[i]]], na.rm=TRUE)
		variance <- var(data[[allvars[i]]], na.rm=TRUE)
		variances <- c(variances, variance)
		if(abs(mean) > 0.1*variance)
			centerwarn <- TRUE
	}	
	if(centerwarn)
		warning('One or more numeric variables has a mean substantially different to 0 - it is highly recommended to centre predictor variables to aid convergence')
	if(length(variances)>1 && any(variances/max(variances, na.rm=TRUE) < 0.01, na.rm=TRUE))
		warning('There is a marked discrepancy in the variance of the numeric predictor variables - it may help convergence to re-scale predictor variables')
	
	# Now set up interactions
	termtypes <- sapply(termvars, function(x){
		stopifnot(length(x)%in%1:2)
		if(length(x)==1)
			return(classes[x])
		if(length(x)==2 && all(c('numeric','factor')%in%classes[x]))
			return('b_int')
		if(length(x)==2 && all(classes[x]=='numeric'))
			return('n_int')
		if(length(x)==2 && all(classes[x]=='factor'))
			return('f_int')
		return(NA)
	})
	
	# First find factor:factor interactions that occur once:
	for(i in which(termtypes=='f_int')){
		# If the two factors don't appear as factor:factor interactions elsewhere then remove the individual factors:
		labels <- termvars[[i]]
		matches <- termtypes=='f_int' & sapply(termvars,function(x) return(any(grepl(labels[1],x)) || any(grepl(labels[2],x))))
		matches[i] <- FALSE
		if(!any(matches)){
			termtypes[i] <- 'f_matrix'
			termtypes[sapply(termvars,length)==1 & sapply(termvars,function(x) return(x[1]))==labels[1]] <- 'dropped_f'
			termtypes[sapply(termvars,length)==1 & sapply(termvars,function(x) return(x[1]))==labels[2]] <- 'dropped_f'
		}
	}
	
	# Then find linear:factor interactions where the linear bit occurs once:
	for(i in which(termtypes=='b_int')){
		# If the linear effect doesn't appear as a factor interaction elsewhere (linear interaction is OK) then drop the main effect:
		labels <- termvars[[i]]
		linvar <- labels[which(classes[labels]=='numeric')]
		matches <- termtypes=='b_int' & sapply(termvars,function(x) return(any(grepl(linvar,x))))
		matches[i] <- FALSE
		if(!any(matches)){
			termtypes[i] <- 'n_matrix'
			termtypes[sapply(termvars,length)==1 & sapply(termvars,function(x) return(x[1]))==linvar] <- 'dropped_n'
		}
	}
	
	varvalues <- vector('list',length=n.chains)
	varnames <- character(0)
	
	if(zifamily){
		extraline <- '\n\tregression_fitted[i] <- regression_positive[i] * non_zero_group[i]'
		respline <- 'regression_positive[i]'
	}else{
		extraline <- ''
		respline <- 'regression_fitted[i]'
	}
	
	extradata <- list(N=nrow(data))
	# Check the response variable:
	if(family=='gaussian'){
		if(!class(data[[response]])%in%c('numeric','integer'))
			stop('The response variable class must be either numeric or integer for the Gaussian family')
		respline <- paste('\t', response, '[i] ~ dnorm(regression_fitted[i], regression_precision)\n\tregression_residual[i] <- ', response, '[i] - regression_fitted[i]\n\tregression_fitted[i] <- ', sep='')
		priorline <- paste('regression_precision ~ ', precision.prior, '\n', sep='')
		varnames <- c(varnames, 'regression_precision')
		signs <- sample(precision.inits, n.chains, replace=TRUE)
		for(c in 1:n.chains)
			varvalues[[c]] <- c(varvalues[[c]], list(signs[c]))		
	}
	if(family=='binomial'){
		if(length(offsets)>0)
			warning('Using an offset() with a logistic regression model is not recommended - did you mean to specify the response variable as a matrix with columns for successes and failures?')
		ok <- FALSE
		if(class(data[[response]])=='matrix'){
			if(ncol(data[[response]])!=2)
				stop('If the response is a matrix, it must have exactly 2 columns')
			data$Binomial_Total <- data[[response]][,1]+data[[response]][,2]
			if(any(is.na(data$Binomial_Total)))
				stop('Missing values are not allowed in the total number of trials')
			if(!all(is.numeric(data$Binomial_Total)) || any(abs(as.integer(data$Binomial_Total)-data$Binomial_Total)>0.001))
				stop('Unexpected non integer value in the numeric response variable (or Binomial_Total variable)')
			data$Binomial_Total <- as.integer(data$Binomial_Total)
			data[[response]] <- data[[response]][,1]
			
			if(!write.data)
				stop('The data must be written to file when supplying a Binomial response as a matrix')
			
		}else{
			if(is.null(data$Binomial_Total))
				data$Binomial_Total <- 1
			
			if(any(is.na(data$Binomial_Total)))
				stop('Missing values are not allowed in the total number of trials')
			if(!all(is.numeric(data$Binomial_Total)) || any(abs(as.integer(data$Binomial_Total)-data$Binomial_Total)>0.001))
				stop('Unexpected non integer value in the supplied Binomial_Total variable')
			data$Binomial_Total <- as.integer(data$Binomial_Total)			
		}
		
		if(all(data$Binomial_Total==1) && zifamily)
			stop('The ZIB model is only available for data with multiple trials - try specifying the data as a matrix, or using Binomial_Total to denote the total number of trials')
		
		if(!write.data && class(data[[response]])%in%c('factor','logical'))
			stop('The data must be written to file when supplying a Binomial response as a factor or logical')
		
    	if(class(data[[response]])=='factor'){
			data[[response]] <- as.numeric(data[[response]])-1
			if(any(data[[response]]>1, na.rm=TRUE))
				warning('Grouping factor levels 2 and above in the response variable')
			ok <- TRUE
		}
		if(class(data[[response]])%in%c('logical','numeric','integer')){
			if(any(abs(as.integer(data[[response]])-data[[response]])>0.001))
				stop('Unexpected non integer value in the numeric response variable (or Binomial_Total variable)')
			data[[response]] <- as.integer(data[[response]])
			if(all(data$Binomial_Total==1) && !all(data[[response]]%in%c(0,1), na.rm=TRUE))
				stop('Unexpected non 0/1 value in the numeric response variable')
			if(any(data$Binomial_Total < data[[response]], na.rm=TRUE))
				stop('Unexpected non 0/1 value in the numeric response variable')
			ok <- TRUE
		}
		priorline <- ''		
		if(!ok)
			stop('Unrecognised response variable format - possibilities are a matrix with columns for successes and failures, or a factor or numeric variable assuming one trial per observation')

		respline <- paste('\t', response, '[i] ~ dbin(regression_prob[i], Binomial_Total[i])\n\tregression_residual[i] <- ', response, '[i] - regression_fitted[i]\n\tregression_fitted[i] <- regression_prob[i] * Binomial_Total[i]\n', if(zifamily) '\tregression_prob[i] <- non_zero_regression[i] * non_zero_group[i]\n\tlogit(non_zero_regression[i]) <- ' else '\tlogit(regression_prob[i]) <- ', sep='')
		extradata <- c(extradata, list(Binomial_Total=data$Binomial_Total))
	}
	if(family=='poisson'){
		if(!class(data[[response]])%in%c('numeric','integer'))
			stop('The response variable class must be either numeric or integer for the Poisson family')
		if(any(data[[response]]<0, na.rm=TRUE) || any(as.integer(data[[response]])!=data[[response]]))
			stop('Only positive integers are allowed in the response variable for the Poisson family')
		respline <- paste('\t', response, '[i] ~ dpois(regression_fitted[i])\n\tregression_residual[i] <- ', response, '[i] - regression_fitted[i]\n\t', if(zifamily) 'regression_fitted[i] <- non_zero_regression[i] * non_zero_group[i]\n\tlog(non_zero_regression[i]) <- ' else 'log(regression_fitted[i]) <- ', sep='')
		priorline <- ''		
	}
	if(family=='nb'){
		if(!class(data[[response]])%in%c('numeric','integer'))
			stop('The response variable class must be either numeric or integer for the Negative Binomial family')
		if(any(data[[response]]<0, na.rm=TRUE) || any(as.integer(data[[response]])!=data[[response]]))
			stop('Only positive integers are allowed in the response variable for the Poisson family')
		respline <- paste('\t', response, '[i] ~ dpois(regression_fitted[i])\n\tregression_residual[i] <- ', response, '[i] - regression_fitted[i]\n\tdispersion[i] ~ dgamma(k, k)\n\tregression_fitted[i] <- regression_mean[i] * dispersion[i]', if(zifamily) ' * non_zero_group[i]', '\n\t# Note: this formulation of a gamma-Poisson is exactly equivalent to a Negative Binomial\n\tlog(regression_mean[i]) <- ', sep='')
		priorline <- paste('k ~ ', precision.prior, '\n\t# Note: the prior for the diserpsion parameter k is quite important for convergence\n\t# [A DuMouchel prior may be better than a Gamma prior]\n', sep='')
		varnames <- c(varnames, 'k')
		signs <- sample(precision.inits, n.chains, replace=TRUE)
		for(c in 1:n.chains)
			varvalues[[c]] <- c(varvalues[[c]], list(signs[c]))		
	}
	if(intercept!=0){
		respline <- paste(respline, 'intercept + ', sep='')
		priorline <- paste(priorline, 'intercept ~ ', effect.prior, '\n', sep='')
		varnames <- c(varnames, 'intercept')
		signs <- sample(effect.inits, n.chains, replace=TRUE)
		for(c in 1:n.chains)
			varvalues[[c]] <- c(varvalues[[c]], list(signs[c]))		
	}
	
	# First numerics, then numeric interactions, then numeric matrices, then factors, then factor matrices, then factor interactions
	for(i in which(termtypes=='numeric')){
		respline <- paste(respline, termvars[[i]][1], '_coefficient * ', termvars[[i]][1], '[i] + ', sep='')
		priorline <- paste(priorline, termvars[[i]][1], '_coefficient ~ ', effect.prior, '\n', sep='')

		varnames <- c(varnames, paste(termvars[[i]][1], '_coefficient', sep=''))
		signs <- sample(effect.inits, n.chains, replace=TRUE)
		for(c in 1:n.chains)
			varvalues[[c]] <- c(varvalues[[c]], list(signs[c]))		
	}
	for(i in which(termtypes=='n_int')){
		respline <- paste(respline, termvars[[i]][1], '_', termvars[[i]][2], '_interaction * ', termvars[[i]][1], '[i] * ', termvars[[i]][2], '[i] + ', sep='')
		priorline <- paste(priorline, termvars[[i]][1], '_', termvars[[i]][2], '_interaction ~ ', effect.prior, '\n', sep='')

		varnames <- c(varnames, paste(termvars[[i]][1], '_', termvars[[i]][2], '_interaction', sep=''))
		signs <- sample(effect.inits, n.chains, replace=TRUE)
		for(c in 1:n.chains)
			varvalues[[c]] <- c(varvalues[[c]], list(signs[c]))		
	}
	for(i in which(termtypes=='n_matrix')){
		linvar <- termvars[[i]][which(classes[termvars[[i]]]=='numeric')]
		fixvar <- termvars[[i]][which(classes[termvars[[i]]]=='factor')]	
		respline <- paste(respline, linvar, '_coefficient_', fixvar, '_level[', fixvar, '[i]] * ', linvar, '[i] + ', sep='')
		factpriors <- rep(paste(' ~ ', effect.prior, sep=''), length(levels(data[[fixvar]])))
		factpriors <- paste(linvar, '_coefficient_', fixvar, '_level[', 1:length(factpriors), ']', factpriors, '    # Factor level "', levels(data[[fixvar]]), '"', sep='')
		priorline <- paste(priorline, paste(factpriors, collapse='\n'), '\n', sep='')

		varnames <- c(varnames, paste(linvar, '_coefficient_', fixvar, '_level', sep=''))
		for(c in 1:n.chains){
			signs <- sample(effect.inits,length(levels(data[[fixvar]])),replace=TRUE)
			varvalues[[c]] <- c(varvalues[[c]], list(signs))
		}
	}
	for(i in which(termtypes=='factor')){
		respline <- paste(respline, termvars[[i]][1], '_effect[', termvars[[i]][1], '[i]] + ', sep='')
		factpriors <- rep(paste(' ~ ', effect.prior, sep=''), length(levels(data[[termvars[[i]][1]]])))
		factpriors[1] <- ' <- 0'
		factpriors <- paste(termvars[[i]][1], '_effect[', 1:length(factpriors), ']', factpriors, '    # Factor level "', levels(data[[termvars[[i]][1]]]), '"', sep='')
		priorline <- paste(priorline, paste(factpriors, collapse='\n'), '\n', sep='')

		varnames <- c(varnames, paste(termvars[[i]][1], '_effect', sep=''))
		for(c in 1:n.chains){
			signs <- sample(effect.inits,length(levels(data[[termvars[[i]][1]]])),replace=TRUE)
			signs[1] <- NA
			varvalues[[c]] <- c(varvalues[[c]], list(signs))
		}
	}
	for(i in which(termtypes=='f_matrix')){
		respline <- paste(respline, termvars[[i]][1], '_', termvars[[i]][2], '_effect[', termvars[[i]][1], '[i],', termvars[[i]][2], '[i]] + ', sep='')
		factindices <- expand.grid(1:length(levels(data[[termvars[[i]][1]]])), 1:length(levels(data[[termvars[[i]][2]]])))
		factpriors <- rep(paste(' ~ ', effect.prior, sep=''), nrow(factindices))
		factpriors[1] <- ' <- 0'
		factpriors <- paste(termvars[[i]][1], '_', termvars[[i]][2], '_effect[', factindices[,1], ',', factindices[,2], '] ', factpriors, '    # Factor level "', levels(data[[termvars[[i]][1]]])[factindices[,1]], '", "', levels(data[[termvars[[i]][2]]])[factindices[,2]], '"', sep='')
		priorline <- paste(priorline, paste(factpriors, collapse='\n'), '\n', sep='')

		varnames <- c(varnames, paste(termvars[[i]][1], '_', termvars[[i]][2], '_effect', sep=''))
		for(c in 1:n.chains){
			signs <- matrix(sample(effect.inits,nrow(factindices),replace=TRUE), nrow=length(levels(data[[termvars[[i]][1]]])))
			signs[1,] <- NA
			signs[,1] <- NA				
			varvalues[[c]] <- c(varvalues[[c]], list(signs))
		}
	}
	for(i in which(termtypes=='f_int')){
		respline <- paste(respline, termvars[[i]][1], '_', termvars[[i]][2], '_interaction[', termvars[[i]][1], '[i],', termvars[[i]][2], '[i]] + ', sep='')
		factindices <- expand.grid(1:length(levels(data[[termvars[[i]][1]]])), 1:length(levels(data[[termvars[[i]][2]]])))
		factpriors <- rep(paste(' ~ ', effect.prior, sep=''), nrow(factindices))
		factpriors[factindices[,1]==1] <- ' <- 0'
		factpriors[factindices[,2]==1] <- ' <- 0'
		factpriors <- paste(termvars[[i]][1], '_', termvars[[i]][2], '_interaction[', factindices[,1], ',', factindices[,2], '] ', factpriors, '    # Factor level "', levels(data[[termvars[[i]][1]]])[factindices[,1]], '", "', levels(data[[termvars[[i]][2]]])[factindices[,2]], '"', sep='')
		priorline <- paste(priorline, paste(factpriors, collapse='\n'), '\n', sep='')

		varnames <- c(varnames, paste(termvars[[i]][1], '_', termvars[[i]][2], '_effect', sep=''))
		signs <- sample(precision.inits, n.chains, replace=TRUE)
		for(c in 1:n.chains){
			signs <- matrix(sample(effect.inits,nrow(factindices),replace=TRUE), nrow=length(levels(data[[termvars[[i]][1]]])))
			signs[1,] <- NA
			signs[,1] <- NA				
			varvalues[[c]] <- c(varvalues[[c]], list(signs))
		}
	}
	# Then random effects:
	for(r in randoms){
		respline <- paste(respline, r, '_randomeffect[', r, '[i]] + ', sep='')
		priorline <- paste(priorline, 'for(', r, '_iterator in 1:', length(levels(data[[r]])), '){\n\t', r, '_randomeffect[', r, '_iterator] ~ dnorm(0, ', r, '_precision)\n}\n', r, '_precision ~ ', precision.prior, '\n', sep='') 
		
		varnames <- c(varnames, paste(r, '_precision', sep=''))
		signs <- sample(precision.inits, n.chains, replace=TRUE)
		for(c in 1:n.chains)
			varvalues[[c]] <- c(varvalues[[c]], list(signs[c]))		
	}
	# And offsets:
	for(o in offsets){
		respline <- paste(respline, o, '[i] + ', sep='')
	}
	
	# Horrible hack to get rid of the trailing +:
	respline <- paste(respline, '_ + _ +\n', sep='')
	respline <- gsub('+ _ + _ +', '', respline, fixed=TRUE)
	
	if(zifamily){
		respline <- paste(respline, '\tnon_zero_group[i] ~ dbern(non_zero_prob[i])\n\tlogit(non_zero_prob[i]) <- -(zero_inflation_intercept)\n\t\t# Note: this line (inside the parentheses) could specify a separate linear regression\n\t\t# To make this the probability of zero-inflation, the - symbol is required!\n', sep='')
		priorline <- paste(priorline, 'zero_inflation_intercept ~ ', effect.prior, '\nnon_zero_propotion <- ilogit(-zero_inflation_intercept)\n', sep='')
		signs <- sample(effect.inits, n.chains, replace=TRUE)
		zistarts <- rep(1, nrow(data))
		varnames <- c(varnames, 'zero_inflation_intercept', 'non_zero_group')
		for(c in 1:n.chains)
			varvalues[[c]] <- c(varvalues[[c]], list(signs[c], zistarts))		
	}
	
	for(c in 1:n.chains)
		names(varvalues[[c]]) <- varnames
	
	# Now create the model file:
	
	if(write.data){
		magicline <- ''
		data <- dump.format(c(data[unique(c(response, allvars, offsets, randoms))], extradata))
	}else{
		magicline <- paste('#data# ', paste(unique(c(response, allvars, offsets, randoms)), collapse=', '), '\n', sep='')
		data <- dump.format(extradata)
	}
	
	if(zifamily){
		modules <- c('dic', 'glm')
	}else{
		modules <- c('glm')
	}
	factories <- ''
	monitor <- c(varnames[!varnames%in%c('non_zero_group','zero_inflation_intercept')], if(zifamily) 'non_zero_propotion' else c('deviance', 'dic'), 'resid.sum.sq')
	
	# Over-write inits with values we have been given:
	for(c in 1:n.chains){
		for(n in names(passedinits[[c]])){
			varvalues[[c]][[n]] <- passedinits[[c]][[n]]
		}
	}
	
	if(write.inits){
		end.state <- sapply(varvalues, dump.format)
	}else{
		end.state <- ''
	}
	model <- paste('### Model template as follows - ensure this is syntactically correct before running the model!\n\nmodel{\n\n# In the BUGS/JAGS language we must use an explicit for loop:\nfor(i in 1:N){\n\t# These lines describe the response distribution and linear model terms:\n', respline, '}\n\n# These lines give the prior distributions for the parameters to be estimated:\n', priorline, 'resid.sum.sq <- sum(regression_residual^2)', magicline, '\n}\n\n# These lines are hooks to be read by runjags (they are ignored by JAGS):', sep='')
	
	rjo <- list(model=model, data=data, end.state=end.state, monitor=monitor, modules=modules, factories=factories, response=response, residual='regression_residual', fitted='regression_fitted')
	class(rjo) <- 'runjags'
	write.jagsfile(runjags.object=rjo, file=file, remove.tags=FALSE)  # Need to keep the data tag in if necessary
	
	# Don't use swcat here deliberately:
	cat('Your model template was created at "', file, '" - it is highly advisable to examine the model syntax to be sure it is as intended\n', sep='')
	if(write.data)
		cat('You can then run the model using run.jags("', file, '")\n', sep='')
	else
		cat('You can then run the model using run.jags("', file, '", data=data) - where "data" is the same data frame specified to the template.jags function\n', sep='')
	
	invisible(file)
	
}


template.JAGS <- template.jags