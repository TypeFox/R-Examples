#
# fit a distribution model to a gene's expression values
#
# - check model fit compared to simple distribution that assumes only measurement error around mean
#
####################################################################################################


fitGMM <- function(exprVals,RDparameters,rejectNull=0.05){

	# RDparameters calculated on non-log-transformed data => take this into account when subtracting 'background' alpha
	Y <- exprVals
	Y[is.na(Y)] <- 0
	
	# compare expression values 'Y' to expected distribution of measurements representing constant expression at median level
	# unfortunately, distribution of error-considering values is difficult to calculate: proceed by drawing data from error distribution

	simulateData <- function(expVal,RDparameters,length){
	
		epsilons <- rnorm(length,mean=0,sd=RDparameters$sd_epsilon)
		epsilons <- (epsilons-mean(epsilons))/sd(epsilons) 
		etas <- rnorm(length,mean=0,sd=RDparameters$sd_eta)
	        etas <- (epsilons-mean(epsilons))/sd(etas)

		y <- 2^expVal
		log((y-(RDparameters$alpha+epsilons))/exp(etas),base=2)
	}

	constantData <- simulateData(median(Y),RDparameters,length=length(Y))
	constantData <- constantData[!is.na(constantData)]
	if(length(constantData)<=2) pval <- 1
	if(length(constantData)>2){
		if(length(constantData)<length(Y)) constantData <- sample(constantData,length(Y),replace=TRUE)
	
		# test that Y,constantData drawn from the same distribution
		# ks test requires no ties
		npoints <- min(length(unique(Y)),length(unique(constantData)))
		kstest <- ks.test(sample(unique(Y),npoints),sample(unique(constantData),npoints))
		pval <- kstest$p.value
		rm(constantData)
	}

	# may be able to decide that genes are varying only due to noise if their variance is less than expected for replicated values with the given mean
	# the 'expected' level of variation could be estimated from resampling values from the error model

	# if test successful to 'rejectNull' p-value threshold, set model to be null:
	if(pval>rejectNull){
		# expression values not significantly differently distributed to simulated measurement error about constant value
		print("not significant differential expression")
		mm <- list(modelName="null",parameters=NULL,G=0)
	}
	else{
		# if test fails, fit Gaussian mixture-model: allow max 5 components and specify non-equal variances
		mm <- Mclust(Y,modelNames=c("V"),G=1:5)
	}
	rm(Y)
	gc()
	mm

}
