`ChiSquaredLRTest` <-
function(
	x,
	model.1.mean.vector, 
	model.1.cov.matrix,
 	model.1.df,
	model.2.mean.vector,
	model.2.cov.matrix,
 	model.2.df
){

	cat(paste("\n", "CHI-SQUARED LR TEST OUTPUT", "\n", sep=""))


	#####################################################################################################################
	# If the data are missing, stop and return an error message
	#####################################################################################################################

	if (missing(x)) {
        	stop("There is no data (x), which must be provided")
    	}



	#####################################################################################################################
	# Subset the data to make a data matrix of only complete cases
	#####################################################################################################################

	x.complete.data <- na.omit(x)



	#####################################################################################################################
	# Extract sample size (N) and the number of measured variables (J)
	# The sample size for complete data (N.complete.data) and the number of cases with missing data (N.missing)

	#####################################################################################################################

	N <- nrow(x)
	J <- ncol(x)

	N.complete.data <- nrow(x.complete.data)
	N.missing = N - N.complete.data




	#####################################################################################################################
	# If the elements of the arguments are missing, use sample values (akin to a saturated model)
	#####################################################################################################################

	if (missing(model.1.mean.vector)) {
        	model.1.mean.vector <- colMeans(x.complete.data)
		cat(paste("\n", "Warning: Sample means will be used for Model 1 mean vector", sep=""))

    	}

	if (missing(model.1.cov.matrix)) {
        	model.1.cov.matrix <- cov(x.complete.data)
		cat(paste("\n", "Warning: Sample covariance matrix will be used for Model 1 covariance matrix", sep=""))
    	}

	if (missing(model.1.df)) {
        	model.1.df <- 0

		cat(paste("\n", "Warning: No df provided for Model 1, analysis will assume 0 df", sep=""))
    	}

	if (missing(model.2.mean.vector)) {
        	model.2.mean.vector <- colMeans(x.complete.data)
		cat(paste("\n", "Warning: Sample means will be used for Model 2 mean vector", sep=""))
    	}

	if (missing(model.2.cov.matrix)) {
        	model.2.cov.matrix <- cov(x.complete.data)
		cat(paste("\n", "Warning: Sample covariance matrix will be used for Model 2 covariance matrix", sep=""))
    	}

	if (missing(model.2.df)) {
        	model.2.df <- 0

		cat(paste("\n", "Warning: No df provided for Model 2, analysis will assume 0 df", sep=""))
    	}

	cat(paste("\n", sep=""))



	#####################################################################################################################
	# If the elements of the arguments are non-conformable, stop and print error
	#####################################################################################################################

	if (NCOL(x.complete.data) != length(model.1.mean.vector)) {
        	stop("Error: Complete data matrix and Model 1 mean vector have non-conforming size")
    	}

	if (NCOL(x.complete.data) != NCOL(model.1.cov.matrix)) {
        	stop("Error: Complete data matrix and Model 1 covariance matrix have non-conforming size")
    	}

	if (NCOL(x.complete.data) != length(model.2.mean.vector)) {
        	stop("Error: Complete data matrix and Model 2 mean vector have non-conforming size")
    	}

	if (NCOL(x.complete.data) != NCOL(model.2.cov.matrix)) {
        	stop("Error: Complete data matrix and Model 2 covariance matrix have non-conforming size")
    	}

	if (NROW(model.1.cov.matrix) != NCOL(model.1.cov.matrix)) {
        	stop("Error: Model 1 covariance matrix must be a square matrix")
    	}

	if (NROW(model.2.cov.matrix) != NCOL(model.2.cov.matrix)) {
        	stop("Error: Model 2 covariance matrix must be a square matrix")
    	}

    	if (length(model.1.mean.vector) != NROW(model.1.cov.matrix)) {
        	stop("Error: Model 1 mean vector and covariance matrix have non-conforming size")
	}

    	if (length(model.2.mean.vector) != NROW(model.2.cov.matrix)) {
        	stop("Error: Model 2 mean vector and covariance matrix have non-conforming size")
	}




	#####################################################################################################################
	# Calculate the T statistic (equation 19 of Levy & Hancock (2007) for comparing models with different BFPDs
	# Begin by initializing arrays
	#####################################################################################################################

	pdf.height.model.1 <- array(NA, N.complete.data)
	pdf.height.model.2 <- array(NA, N.complete.data)



	#####################################################################################################################
	# For each subject with complete data calculate the height of the pdf under each model
	#####################################################################################################################

	for(i in 1:N.complete.data){

		pdf.height.model.1[i] <- dmvnorm(x=x.complete.data[i,], mean=	model.1.mean.vector, sigma=model.1.cov.matrix, log=FALSE) 
		pdf.height.model.2[i] <- dmvnorm(x=x.complete.data[i,], mean=	model.2.mean.vector, sigma=model.2.cov.matrix, log=FALSE) 

	
	}



	#####################################################################################################################
	# Calculate the individual level LR
	# Calculate the LR 
	# Calculate the test statistic T and the 2-tailed p-value
	#####################################################################################################################

	log.likelihood.individ <- log(pdf.height.model.1)-log(pdf.height.model.2)
	LR <- sum(log.likelihood.individ)

	chi.squared.diff <- -2*LR
	df.diff <- model.1.df-model.2.df

	p = pchisq(q=chi.squared.diff, df=df.diff, lower.tail=FALSE)



	#####################################################################################################################
	# Print the results to the screen
	#####################################################################################################################

	cat(paste("\n", sep=""))
	cat(paste("RESULTS OF CHI-SQUARED LR TEST", "\n", "\n", sep=""))
	cat(paste("Total sample size = ", N, "\n", sep=""))
	cat(paste("Sample size with complete data used in the analysis = ", N.complete.data, "\n", sep=""))
	cat(paste("LR = ", round(LR, 3), "\n", sep=""))
	cat(paste("Chi-Squared Statistic = ", round(chi.squared.diff, 3), "\n", sep=""))
	cat(paste("df for the test = ", df.diff, "\n", sep=""))
	cat(paste("p-value = ", round(p,3), "\n", sep=""))
	cat(paste("\n", sep=""))



	#####################################################################################################################
	# Return the selected results
	#####################################################################################################################

	list(N = N, N.complete.data = N.complete.data, LR=LR, chi.sq.stat=chi.squared.diff, df.test=df.diff, p = p)


}

