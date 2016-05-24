##################################################################################
# Set Method: classifyRasclass
##################################################################################
classifyRasclass <- function(rasclassObj, splitfraction = 1, method = 'logit', ...){}
setMethod('classifyRasclass', signature(rasclassObj = 'rasclass'),

function(rasclassObj, splitfraction = 1, method = 'logit', ...){
	
	# Control for Forumla
	if(is.null(rasclassObj@formula)){
		stop('Please build or specify formula in rasclass object')
	}
	
	# Check consistency
	if(!checkRasclass(rasclassObj)){
		stop('Data object is not consistent, check data and try again')
	}
	
	# Store call
	rasclassObj@call <- match.call()
	
	##################################################################################
	# Create index to randomly split data in two parts
	
	if((splitfraction < 1) & (splitfraction > 0)){
		
		samplesize <- sum(!is.na(rasclassObj@data[, rasclassObj@samplename]))
		separator   <- sample(c(
			rep(TRUE, round(samplesize*splitfraction)),
			rep(FALSE, samplesize - round(samplesize*splitfraction))
			))

		training <- !is.na(rasclassObj@data[, rasclassObj@samplename])
		
		training[training] <- separator

		rasclassObj@training <- as.logical(training)
		
	} else if(splitfraction != 1){
		stop('Splitfraction not in interval (0, 1]')
	}	
	
	##################################################################################
	# Classify the dataset using specified classification algorithm

	# Make sure the sample is stored as factor
	if(class(rasclassObj@data[, rasclassObj@samplename]) != 'factor'){
		rasclassObj@data[, rasclassObj@samplename] <- factor(rasclassObj@data[, rasclassObj@samplename])
	}

	# Prepare predictor data list (it is more efficient to predict in pieces)
	if(method != 'maximumLikelihood'){
		splitter <- factor(rep(1:(1 + nrow(rasclassObj@data)/1e4), length.out = nrow(rasclassObj@data), each = 1e4))
		predictorVars <- split(rasclassObj@data[, names(rasclassObj@data) != rasclassObj@samplename], splitter)			
	}
	
	# Switch by Method type
	######### Multinomial Logistic Regression #########
	if(method == 'logit'){
		cat('Classification using Multinomial Logistic Regression\n')

		cat('classifying...\n')
		if(splitfraction != 1){
			rasclassObj@logit <- multinom(rasclassObj@formula, rasclassObj@data[training, ], ...)
		} else {
			rasclassObj@logit <- multinom(rasclassObj@formula, na.omit(rasclassObj@data), ...)
		}
		cat('predicting...\n')

		predicted <- lapply(predictorVars, function(x){
			cat('|')
			as.character(predict(rasclassObj@logit, x))
			})
		predicted <- as.numeric(do.call(c, predicted))
	}
	######### randomForest #########
	else if(method == 'randomForest'){
		cat('Classification using Random Forest\n')

		cat('classifying...\n')
		if(splitfraction != 1){
			rasclassObj@randomForest <- randomForest(as.formula(rasclassObj@formula), rasclassObj@data[training, ], ...)
		} else {
			rasclassObj@randomForest <- randomForest(as.formula(rasclassObj@formula), na.omit(rasclassObj@data), ...)
		}
		cat('predicting...\n')
				predicted <- lapply(predictorVars, function(x){
			cat('|')
			as.character(predict(rasclassObj@randomForest, x))
			})
		predicted <- as.numeric(do.call(c, predicted))
	}
	######### Support Vector Machines #########
	else if(method == 'supportVector'){
		cat('Classification using Support Vector Machines\n')

		cat('classifying...\n')
		if(splitfraction != 1){
			rasclassObj@supportVector <- svm(as.formula(rasclassObj@formula), rasclassObj@data[training, ], ...)
		} else {
			rasclassObj@supportVector <- svm(as.formula(rasclassObj@formula), na.omit(rasclassObj@data), ...)
		}
		cat('predicting...\n')
		predicted <- lapply(predictorVars, function(x){
			cat('|')
			as.character(predict(rasclassObj@supportVector, x, decision.values = TRUE))
			})
		predicted <- as.numeric(do.call(c, predicted))
	}
	######### Neural Network #########
	else if(method == 'neuralNetwork'){
		cat('Classification using Neural Network\n')
		
		# Create required Lists and varaibles
		varlist <- strsplit(rasclassObj@formula,' ~ ', TRUE)[[1]]
		varlist <- unlist(strsplit(as.character(varlist[2]),' + ', TRUE))

		cat('classifying...\n')
		if(splitfraction != 1){
			classLabels <- matrix(0,
				nrow = length(rasclassObj@data[training, rasclassObj@samplename]),
				ncol = nlevels(rasclassObj@data[training, rasclassObj@samplename]),
				dimnames = list(NULL, levels(rasclassObj@data[training, rasclassObj@samplename])))
			sampleIndex <-as.integer(rasclassObj@data[training, rasclassObj@samplename])
			for(row in 1:nrow(classLabels)){
				classLabels[row, sampleIndex[row]] <- 1
			}

			normedData <- normalizeData(rasclassObj@data[training, varlist])
			rasclassObj@neuralNetwork <- mlp(x = normedData, y = classLabels)
		} else {
			classLabels <- matrix(0,
				nrow = length(na.omit(rasclassObj@data)[, rasclassObj@samplename]),
				ncol = nlevels(na.omit(rasclassObj@data)[, rasclassObj@samplename]),
				dimnames = list(NULL, levels(na.omit(rasclassObj@data)[, rasclassObj@samplename])))
			sampleIndex <-as.integer(na.omit(rasclassObj@data)[, rasclassObj@samplename])
			for(row in 1:nrow(classLabels)){
				classLabels[row, sampleIndex[row]] <- 1
			}
			
			normedData <- normalizeData(na.omit(rasclassObj@data)[varlist])			
			rasclassObj@neuralNetwork <- mlp(x = normedData, y = classLabels)
		}
		cat('predicting...\n')
		predictorVars <- lapply(predictorVars, function(x) x[varlist])
		predicted <- lapply(predictorVars, function(x){
			cat('|')
			nn.pred <- predict(rasclassObj@neuralNetwork, normalizeData(x, type = getNormParameters(normedData)))
			nn.pred <- as.numeric(colnames(classLabels)[apply(nn.pred, 1, which.max)])
			nn.pred
		})
		predicted <- do.call(c, predicted)
	}
	######### MLC #########
	else if(method == 'maximumLikelihood'){
		cat('Classification using Maximum Likelihood Classifier\n')
		mlc.out <- rasclassMlc(rasclassObj)
		rasclassObj@maximumLikelihood <- mlc.out[[1]]
		predicted <- mlc.out[[2]]
	}
	######### Falsely specified classifier #########
	else{
		stop(paste('Method', method, 'is not defined in function classifyRasclass. Allowed methods are: \'maximumLikelihood\', \'logit\', \'neuralNetwork\', \'randomForest\' and \'supportVector\'.'))
	}

	# Remove predictor list
	if(method != 'maximumLikelihood'){
		rm(predictorVars)
		gc()
	}
	
	##################################################################################
	# Store predicted Grid using skeleton
	cat('\nStoring predicted grid...\n')
	rasclassObj@predictedGrid@ncols     <- rasclassObj@gridSkeleton@ncols
	rasclassObj@predictedGrid@nrows     <- rasclassObj@gridSkeleton@nrows
	rasclassObj@predictedGrid@xllcorner <- rasclassObj@gridSkeleton@xllcorner
	rasclassObj@predictedGrid@yllcorner <- rasclassObj@gridSkeleton@yllcorner
	rasclassObj@predictedGrid@cellsize  <- rasclassObj@gridSkeleton@cellsize
	rasclassObj@predictedGrid@NAvalue   <- rasclassObj@gridSkeleton@NAvalue
	rasclassObj@predictedGrid@grid      <- as.numeric(rep(NA, length(rasclassObj@gridSkeleton@grid)))
	
	rasclassObj@predictedGrid@grid[as.logical(rasclassObj@gridSkeleton@grid)] <- predicted

	##################################################################################
	# Calculate accuracy of predicted grid
	cat('Calculating accuracy Matrix...\n')
	# Confusion matrix
	if(splitfraction != 1){
		accuracytable <- as.matrix(table(
			rasclassObj@data[, rasclassObj@samplename][!training],
			factor(predicted[!training], levels = levels(rasclassObj@data[, rasclassObj@samplename][!training]))))

	} else {
		accuracytable <- as.matrix(table(
			rasclassObj@data[, rasclassObj@samplename],
			factor(predicted, levels = levels(rasclassObj@data[, rasclassObj@samplename]))))
	}	
	
	# Adapt names of matrix
	rownames(accuracytable) <- paste('Sample', rownames(accuracytable))
	colnames(accuracytable) <- paste('Predicted', colnames(accuracytable))
	
	# Overall accuracy
	rasclassObj@overallAccuracy <- sum(diag(accuracytable))/sum(accuracytable)
	
	# Kappa coefficient	
	thissum = sum(rowSums(accuracytable)*colSums(accuracytable))
	n <- sum(rowSums(accuracytable))
	diag <- sum(diag(as.matrix(accuracytable)))
	rasclassObj@kappa <- (n*diag - thissum)/(n^2 - thissum)
	
	# User's and Producer's accuracy
	produceraccuracy <- as.matrix(accuracytable/rowSums(accuracytable))
	useraccuracy     <- as.matrix(accuracytable/colSums(accuracytable))
	
	accuracytable <- cbind(accuracytable, diag(produceraccuracy))
	colnames(accuracytable)[ncol(accuracytable)] <- 'Producer Acc'

	accuracytable <- rbind(accuracytable,append(diag(useraccuracy),c(NA)))
	rownames(accuracytable)[nrow(accuracytable)] <- 'User Acc'
	
	rasclassObj@accuracyMatrix <- accuracytable

	rasclassObj
}
)