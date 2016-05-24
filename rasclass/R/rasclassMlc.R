##################################################################################
# Set Method: rasclassMLC
##################################################################################
rasclassMlc <- function(rasclassObj){}
setMethod('rasclassMlc', signature(rasclassObj = 'rasclass'), definition = 

function(rasclassObj){

	# Create required Lists and varaibles
	varlist    <- strsplit(as.character(rasclassObj@formula),' ~ ', TRUE)[[1]]
	samplename <- varlist[1]
	varlist    <- unlist(strsplit(as.character(varlist[2]),' + ', TRUE))
	
	# Check consistency of input
	if((sum(sapply(varlist, function(x) !is.element(x, names(rasclassObj@data)))) != 0) | (!is.element(samplename, names(rasclassObj@data)))){
		stop('Formula not consistent with data frame names')
	}
	
	# Calculate parameters of the multivariate normal distribution
	cat('classifying...\n')
	coefs <- list()
	if(is.null(rasclassObj@training)){
		byClass <- split(rasclassObj@data[, varlist], rasclassObj@data[, samplename])
		samplesize <- length(na.omit(rasclassObj@data[, samplename]))
	} else {
		byClass <- split(rasclassObj@data[rasclassObj@training, varlist], rasclassObj@data[rasclassObj@training, samplename])
		samplesize <- sum(rasclassObj@training)
	}
	
	classes <- as.numeric(names(byClass))
	for(cat in names(byClass)){
		frame <- byClass[[cat]]
		
		# Calculate parameters of the multivariate normal distribution
		prior <- log(nrow(frame)/samplesize)
		meanVector  <- colMeans(frame)
		classCov <- cov(frame)
		if(sum(diag(classCov) == 0) != 0){
			failnames <- names(diag(classCov))[diag(classCov) == 0]
			warning('No variation of variable(s) "', paste(failnames, collapse = ', '), '" within class ', cat,'\n. Ignoring variable for prediction in this class.')
			diag(classCov)[diag(classCov) == 0] <- 1
		}
		determinant <- log(det(classCov))
		inverseCov  <- solve(classCov)
		coefs[[cat]] <- list(prior, determinant, meanVector, inverseCov)
	}
	rm(byClass)
	gc()

	# Predict grid values based on probabilities
	cat('predicting...\n')
	dataVars <- as.matrix(rasclassObj@data[, varlist])
	predicted <- rep(NA, nrow(dataVars))
	probs <- rep(NA, length(classes))

	for(i in 1:nrow(dataVars)){
		for(j in 1:length(probs)){
			delta <- dataVars[i,] - coefs[[j]][[3]]
			probs[j] <- coefs[[j]][[1]] - coefs[[j]][[2]] - t(delta)%*%coefs[[j]][[4]]%*%delta
		}
		predicted[i] <- classes[probs==max(probs)]

		if(i%%10000 == 1) cat('|')
	}
	
	#Return predicted vector
	out <- list()
	out[[1]] <- coefs
	out[[2]] <- predicted
	out
}
)
