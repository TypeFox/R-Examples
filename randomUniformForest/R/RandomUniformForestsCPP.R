# <OWNER> = Saip CISS
# <ORGANIZATION> = QueensBridge Quantitative
# <YEAR> = 2014

# LICENSE 
# BSD 3-CLAUSE LICENSE

# Copyright (c) 2014, Saip CISS
# All rights reserved.

# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

# 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# END OF LICENSE 

randomUniformForest <- function(...) UseMethod("randomUniformForest")
	
importance <- function(object, ...) UseMethod("importance")

unsupervised <- function(object,...) UseMethod("unsupervised")

randomUniformForest.formula <- function(formula, data = NULL, subset = NULL, ...)
{
	if (is.null(data))stop ("Please provide data.\n")	
	if (!is.null(subset)) data = data[subset,]	
	data <- fillVariablesNames(data)	
	mf <- model.frame(formula = formula, data = as.data.frame(data))
	x <- model.matrix(attr(mf, "terms"), data = mf)[,-1]
	y <- model.response(mf)
	names(y) = NULL
	RUFObject <- randomUniformForest.default(x, Y = y, ...)
	RUFObject$call <- match.call()
	RUFObject$formula <- formula
	class(RUFObject) <- c("randomUniformForest.formula", "randomUniformForest")
	RUFObject
}

print.randomUniformForest <- function(x,...)
{
	object <- x
	cat("Call:\n")
	print(object$call)
	cat("\n")
	cat("Type of random uniform forest: ") 
	if (!is.null(object$unsupervised))
	{	cat("Unsupervised learning\n") }
	else
	{
		if (object$forest$regression) 	{  cat("Regression\n")	}
		else { 	cat("Classification\n") }
	}
	cat("\n")
	print(object$forestParams)
	cat("\n")
	if (!is.null(object$forest$OOB))
	{ 
		cat("Out-of-bag (OOB) evaluation")
		if (!object$forest$regression)
		{
			if (!is.numeric(object$forest$pred.error))
			{	cat("\n", object$forest$pred.error, "\n", "Options used seem have to be reconsidered.","\n") }
			else
			{
				OOBErrorRate = mean(100*object$forest$pred.error)
				cat("\nOOB estimate of error rate: ", round(OOBErrorRate, 2),"%\n", sep ="")
				cat("OOB error rate bound (with 1% deviation): ",round(OOBErrorRate + OOBErrorRate*(1 - estimatePredictionAccuracy(floor(length(object$forest$OOB.predicts)*0.368))), 2),"%\n", sep = "")
				cat("\nOOB confusion matrix:\n")
				colnames(object$forest$OOB)[1:length(object$classes)] = object$classes
				rownames(object$forest$OOB)[1:length(object$classes)] = object$classes
				print(round(object$forest$OOB,4))			
				if ((length(object$classes) == 2) & (rownames(object$forestParams)[1] != "reduceDimension"))
				{	
					cat("\nOOB estimate of AUC: ", round(pROC::auc(as.numeric(object$y), as.numeric(object$forest$OOB.predicts))[[1]], 4), sep = "")
					cat("\nOOB estimate of AUPR: ", round(myAUC(as.numeric(object$forest$OOB.predicts), 
					as.numeric(object$y), falseDiscoveryRate = TRUE)$auc, 4),sep = "")
					cat("\nOOB estimate of F1-score: ", round(fScore(object$forest$OOB), 4),sep = "")
				}
				cat("\nOOB (adjusted) estimate of geometric mean:", round(gMean(object$forest$OOB),4),"\n")
				if (nrow(object$forest$OOB) > 2)
				{ cat("OOB (adjusted) estimate of geometric mean for precision:", round(gMean(object$forest$OOB, precision = TRUE),4),"\n") }
				if (!is.null(object$forest$OOB.strengthCorr))
				{
					cat("\nBreiman's bounds")
					cat("\nExpected prediction error (under approximatively balanced classes): ", round(mean(100*rmNA(object$forest$OOB.strengthCorr$PE)),2),"%\n", sep ="")
					cat("Upper bound: ", round(mean(100*rmNA(object$forest$OOB.strengthCorr$std.strength^2/object$forest$OOB.strengthCorr$strength^2)),2),"%\n", sep ="")
					cat("Average correlation between trees:", round(mean(round(rmNA(object$forest$OOB.strengthCorr$avg.corr),4)),4),"\n")
					cat("Strength (margin):", round(mean(round(rmNA(object$forest$OOB.strengthCorr$strength),4)),4),"\n")
					cat("Standard deviation of strength:", round(mean(round(rmNA(object$forest$OOB.strengthCorr$std.strength),4)),4),"\n")
				}
			}
		}
		else
		{
			cat("\nMean of squared residuals:", round(mean(object$forest$pred.error),8), "\n")
			cat("Mean squared error bound (experimental):", round(object$forest$pred.error + object$forest$pred.error*(1- estimatePredictionAccuracy(floor(length(object$forest$OOB.predicts)*
			(1 - as.numeric(as.vector(object$forestParams["subsamplerate",1])))))), 6),"\n")
			if (length(object$forest$percent.varExplained) > 1)
			{	
				varExplained = object$forest$percent.varExplained
				for (i in 1:length(varExplained))
				{	varExplained[i] = paste(object$forest$percent.varExplained[i], "%", sep="")	}
				cat("Variance explained:", varExplained, "\n")	
			}
			else { cat("Variance explained: ", object$forest$percent.varExplained, "%\n", sep = "")	} #percent.varExplained
			cat("\nOOB residuals:\n")
			Residuals = summary(rmNA(object$forest$OOB.predicts - object$y))
			names(Residuals) = c("Min", "1Q", "Median", "Mean", "3Q", "Max")
			print(Residuals)
			cat("Mean of absolute residuals:", sum(abs(rmNA(object$forest$OOB.predicts - object$y)))/length(rmNA(object$y)),"\n")			
			if (!is.null(object$forest$OOB.strengthCorr))
			{
				cat("\nBreiman's bounds")
				cat("\nTheoretical prediction error:", round(rmNA(object$forest$OOB.strengthCorr$PE.forest),6) ,"\n")
				cat("Upper bound:", round(rmNA(object$forest$OOB.strengthCorr$PE.max),6),"\n")
				cat("Mean prediction error of a tree:", round(rmNA(object$forest$OOB.strengthCorr$PE.tree),6),"\n")
				cat("Average correlation between trees residuals:", round(rmNA(object$forest$OOB.strengthCorr$mean.corr),4),"\n")
				residuals_hat = vector(length = ncol(object$forest$OOB.votes))
				residuals_hat <- apply(object$forest$OOB.votes, 2, function(Z) mean(rmInf(Z) - mean(rmNA(object$forest$OOB.predicts))))		
				cat("Expected squared bias (experimental):", round(mean(rmNA(residuals_hat))^2,6),"\n")
			}
		}		
	}  
    if (!is.null(object$errorObject))
	{	
	    cat("\nTest set")
		if (!object$forest$regression)
		{
			cat("\nError rate: ")
			cat(round(100*object$errorObject$error,2), "%\n", sep="")
			cat("\nConfusion matrix:\n")
			print(round(object$errorObject$confusion,4))			
			if (!is.null(object$errorObject$AUC)) 
			{ 
				cat("\nArea Under ROC Curve:", round(object$errorObject$AUC,4))
				cat("\nArea Under Precision-Recall Curve:", round(object$errorObject$AUPR,4))
				cat("\nF1 score:", round(fScore(object$errorObject$confusion),4)) 	
			}
			cat("\nGeometric mean:", round(gMean(object$errorObject$confusion),4),"\n")
			if (nrow(object$errorObject$confusion) > 2)
			{ cat("Geometric mean of the precision:", round(gMean(object$errorObject$confusion, precision = TRUE),4),"\n") }
		}
		else
		{	
			cat("\nMean of squared residuals: ", round(object$errorObject$error, 6), "\n", sep="")			
			cat("Variance explained: ", object$errorObject$percent.varExplained, "%\n\n", sep = "")		
			Residuals <- object$errorObject$Residuals
			names(Residuals) = c("Min", "1Q", "Median", "Mean", "3Q", "Max")
			cat("Residuals:\n")
			print(Residuals)
			cat("Mean of absolute residuals: ", round(rmNA(object$errorObject$meanAbsResiduals), 6), "\n", sep="")		
		} 		
	}
}

summary.randomUniformForest <- function(object, maxVar = 30, border = NA,...)
{	
	object <- filter.object(object)	
	if (!is.null(object$forest$variableImportance))
	{	
		par(las=1) 
		maxChar = floor(2 + max(nchar(object$variablesNames))/2)		
		par(mar=c(5, maxChar + 1,4,2))		
		varImportance1 = varImportance = object$forest$variableImportance		
		if (!object$forest$regression)
		{
			varImportance[,"class"] = object$classes[as.numeric(varImportance[,"class"])]
			varImportance1[,"class"] = varImportance[,"class"]
		}		
		nVar = nrow(varImportance)
		if (nVar > maxVar) { varImportance = varImportance[1:maxVar,] }	
		barplot( varImportance[nrow(varImportance):1,"percent.importance"], horiz = TRUE, 
		col = sort(heat.colors(nrow(varImportance)), decreasing = TRUE), 
		names.arg = varImportance[nrow(varImportance):1,"variables"], border = border,
		xlab = "Relative Influence (%)", 
		main = 	if (!object$forest$regression) 
				{ "Variable Importance based on Information Gain" }
				else 
				{ 
					if ( length(grep(as.character(object$forestParams[nrow(object$forestParams),1]), "absolute")) == 1 ) 
					{ "Variable Importance based on L1 distance" }
					else
					{ "Variable Importance based on L2 distance" }
				})
		abline(v = 100/nVar, col ='grey')
		cat("\nGlobal Variable importance:\n")
		if (!object$forest$regression)
		{
			cat("Note: most predictive features are ordered by 'score' and plotted. Most discriminant ones\nshould also be taken into account by looking 'class' and 'class.frequency'.\n\n")
		}
		print(varImportance1)
		cat("\n")
		cat("Average tree size (number of nodes) summary: ", "\n")
		nodesView = unlist(lapply(object$forest$object,nrow))
		print(floor(summary(nodesView)))
		cat("\n")
		cat("Average Leaf nodes (number of terminal nodes) summary: ", "\n")
		terminalNodesView = unlist(lapply(object$forest$object, function(Z) length(which(Z[,"status"] == -1))))
		print(floor(summary(terminalNodesView)))
		cat("\n")
		cat("Leaf nodes size (number of observations per leaf node) summary: ", "\n")
		print(summary(unlist(lapply(object$forest$object, function(Z) Z[which(Z[,"prediction"] != 0), "nodes"]) )))
		cat("\n")
		cat("Average tree depth :", round(log(mean(nodesView))/log(2),0), "\n")
		cat("\n")
		cat("Theoretical (balanced) tree depth :", round(log(length(object$y), base = 2),0), "\n")
		cat("\n")		
	}
	else
	{ 
		cat("Average tree size (number of nodes) summary: ", "\n")
		nodesView = unlist(lapply(object$forest$object,nrow))
		print(floor(summary(nodesView)))
		cat("\n")
		cat("Average Leaf nodes (number of terminal nodes) summary: ", "\n")
		terminalNodesView = unlist(lapply(object$forest$object, function(Z) length(which(Z[,"status"] == -1))))
		print(floor(summary(terminalNodesView)))
		cat("\n")
		cat("Leaf nodes size (number of observations per leaf node) summary: ", "\n")
		print(summary(unlist(lapply(object$forest$object, function(Z) Z[which(Z[,"prediction"] != 0), "nodes"]) )))
		cat("\n")
		cat("Average tree depth :", round(log(mean(nodesView), base = 2),0), "\n")
		cat("\n")
		cat("Theoretical (balanced) tree depth :", round(log(length(object$y), base = 2),0), "\n")
		cat("\n")		
	}
}

plot.randomUniformForest <- function(x, threads = "auto", ...) 
{
	object <- x
	if (!is.null(object$forest$OOB.votes))
	{
		Ytrain = object$y
		OOBMonitoring = object$forest$OOB.votes
		
		ff = L2Dist
		ffComment = "OOB mean squared error"
		if ( length(grep(as.character(object$forestParams[nrow(object$forestParams),1]), "absolute")) == 1 ) 
		{ 
			ff = L1Dist	
			ffComment = "OOB mean absolute error"
		}
				
		ZZ <- monitorOOBError(OOBMonitoring, Ytrain, regression = object$forest$regression, threads = threads, f = ff)
		
		if (object$forest$regression)
		{	plot(ZZ, type = 'l', lty=2, xlab = "Trees", ylab = ffComment, ...) 	}
		else
		{	    
			plot(apply(ZZ[,1:3],1, min), type='l', lty=2, col = "green", xlab = "Trees", ylab ="OOB error", ...)
			points(apply(ZZ[,1:3],1, mean), type='l', lty=3)
			points(apply(ZZ[,1:3],1, max), type='l', lty=3, col='red')
		}
		grid()
	}
	else
	{ print("no OOB data to plot")	}
}

getTree.randomUniformForest <- function(object, whichTree, labelVar = TRUE)
{
	if (labelVar)
	{
		Tree = data.frame(object$forest$object[[whichTree]])
		idx = which(Tree[, "split.var"] != 0)
		Tree[idx, "split.var"] = object$variablesNames[Tree[idx, "split.var"]]
		
		return(Tree)
	}
	else
	{  return(object$forest$object[[whichTree]])  }
}
		
predict.randomUniformForest <- function(object, X, 
	type = c("response", "prob", "votes", "confInt", "ranking", "quantile", "truemajority", "all"),
	classcutoff = c(0,0), 
	conf = 0.95,
	whichQuantile = NULL,
	rankingIDs = NULL,
	threads = "auto", 
	parallelpackage = "doParallel", ...) rUniformForestPredict(object, X, type = type, classcutoff = classcutoff, 
									conf = conf,
									whichQuantile = whichQuantile,
									rankingIDs = rankingIDs,
									threads = threads, 
									parallelpackage = parallelpackage, ...)

residualsRandomUniformForest <- function(object, Y = NULL) 
{
	object = filter.object(object)
	if (is.null(Y))
	{
		if (is.null(object$forest$OOB.predicts))
		{	stop("please enable OOB option when computing a random uniform forest") }
		else
		{	
			print("OOB residuals:")
			cat("\n")
			return(object$y - object$forest$OOB.predicts) 
		}
	}
	else
	{	
		if (is.numeric(object))
		{ return(object - Y) }
		else
		{ 
			if (!is.null(object$predictionObject$majority.vote))
			{ return(object$predictionObject$majority.vote - Y) }
			else
			{	stop("Please provide model responses to compute residuals")	}
		}
	}	
}

genericOutput <- function(xtest, ytest, paramsObject, RUF.model, ytrain = NULL, classcutoff = c(0,0))
{
	classes = NULL
	if (!is.null(ytest))
	{
		if (is.factor(ytest)) { YNames = classes = levels(ytest)   }
	}
	else
	{
		if (!is.null(ytrain))
		{
			if (!RUF.model$regression) 
			{ 
				ytrain = as.factor(ytrain) 
				YNames = classes = levels(ytrain)  
			}
		}
		else
		{
			if (!RUF.model$regression)  { YNames = classes = sort(unique(RUF.model$object[[2]][,6])[-1]) }
		}
	}		
	if (!is.null(RUF.model$OOB))
	{
		if (!RUF.model$regression & is.factor(ytest) & !is.character(RUF.model$OOB)) 
		{ row.names(RUF.model$OOB) = colnames(RUF.model$OOB)[-(length(YNames)+1)] = YNames }
	}
	if (!is.null(xtest)) 
	{ 
		classwtString = as.vector(paramsObject[which(row.names(paramsObject) == "classwt"),1])	
		classwt = FALSE
		if (is.na(as.logical(classwtString))) {  classwt = TRUE }
		
		if (!RUF.model$regression & (as.numeric(classcutoff[2]) != 0)) 
		{  
			classcutoff = c(which(classes == as.character(classcutoff[1])), as.numeric( classcutoff[2]))
		}				
		if (classcutoff[2] != 0 ) {  classcutoff[2] = 0.5/classcutoff[2] }		
		RUF.predicts <- randomUniformForestCore.predict(RUF.model, xtest, pr.classwt = classwt, 
		pr.imbalance = classcutoff)		 
		if (!is.null(ytest)) 
		{ 
			majorityVote = RUF.predicts$majority.vote
			if (!RUF.model$regression)
			{
				majorityVote = as.factor(majorityVote) 
				levels(majorityVote) = classes[as.numeric(levels(majorityVote))]				
			}
			errorObject <- someErrorType(majorityVote, ytest, regression = RUF.model$regression)
			if (!RUF.model$regression & is.factor(ytest)) 
			{ row.names(errorObject$confusion) = colnames(errorObject$confusion)[-(length(YNames)+1)] = YNames }			
			RUFObject = list(forest = RUF.model, predictionObject = RUF.predicts, errorObject = errorObject, forestParams = paramsObject, classes = classes) 
		}
		else
		{	
			RUFObject = list(forest = RUF.model, predictionObject = RUF.predicts, forestParams = paramsObject,  
			classes = classes) 
		}		
	}
	else
	{   RUFObject = list(forest = RUF.model, forestParams = paramsObject, classes = classes)	}	
	RUFObject
}

randomUniformForest.default <- function(X, Y = NULL, xtest = NULL, ytest = NULL, ntree = 100, 
	mtry = ifelse(bagging,ncol(X),floor(4/3*ncol(X))), 
	nodesize = 1,
	maxnodes = Inf,
    depth = Inf,
    depthcontrol = NULL,	    
	regression = ifelse(is.factor(Y), FALSE, TRUE),
	replace = ifelse(regression,FALSE,TRUE),
	OOB = TRUE,
	BreimanBounds = ifelse(OOB, TRUE, FALSE),
	subsamplerate = ifelse(regression,0.7,1),
    importance = TRUE,	
	bagging = FALSE,
	unsupervised = FALSE, 
	unsupervisedMethod = c("uniform univariate sampling", "uniform multivariate sampling", "with bootstrap"),
	classwt = NULL,
	oversampling = 0,
	targetclass = -1,
	outputperturbationsampling = FALSE,
    rebalancedsampling = FALSE,
	featureselectionrule = c("entropy", "gini", "random", "L2", "L1"),	
	randomcombination = 0,
	randomfeature = FALSE,
	categoricalvariablesidx = NULL,
	na.action = c("fastImpute", "accurateImpute", "omit"),
	logX = FALSE,
	classcutoff = c(0,0),
	subset = NULL,
	usesubtrees = FALSE,
	threads = "auto",
	parallelpackage = "doParallel", 
	...)
{
	{
		if (threads != 1)
		{
			if (sample(9,1) == 9) { rm.tempdir() }
		}
		if (!is.null(subset))
		{
			X = X[subset,]
			Y = Y[subset]
		}
		if (exists("categorical", inherits = FALSE)) { categoricalvariablesidx = categorical }
		else { categorical = NULL  }
		if (is.null(Y) | (unsupervised == TRUE)) 
		{ 
			cat("Enter in unsupervised learning.\n")
			if (unsupervisedMethod[1] == "with bootstrap")
			{ 
				cat("'with bootstrap' is only needed as the second argument of 'unsupervisedMethod' option. Default option will be computed.\n")
				unsupervisedMethod[1] = "uniform univariate sampling"
			}
			XY <- unsupervised2supervised(X, method = unsupervisedMethod[1], bootstrap = if (length(unsupervisedMethod) > 1) { TRUE } else {FALSE })
			X = XY$X
			Y = as.factor(XY$Y)
			unsupervised = TRUE
		}
		if (ntree < 2) { stop("Please use at least 2 trees for computing forest.\n") }			
		if ( (subsamplerate == 1) &  (replace == FALSE))  { OOB = FALSE }
		if (subsamplerate > 1) { replace = TRUE }
		if (depth < 3) { stop("Stumps are not allowed. Minimal depth is 3, leading, at most, to 8 leaf nodes.\n") }
		if (!is.null(classwt) & (!is.factor(Y) | regression)) 
		{ 
			cat("Class reweighing is not allowed for regression. Resetting to default values.\n")
			classwt = NULL
		}
		if (targetclass == 0) { stop("'targetclass' must take a strictly positive value") }
		if ( (BreimanBounds & (length(unique(Y)) > 2) & (ntree > 500)) | (BreimanBounds & (ntree > 500)) ) 
		{ cat("Note: Breiman's bounds (especially for multi-class problems) are computationally intensive.\n") }
		if (maxnodes < 6) { stop("Maximal number of nodes must be above 5.\n") }
		X <- fillVariablesNames(X)        
		getFactors <- which.is.factor(X, count = TRUE)
		if ( (sum(getFactors) > 0) & is.null(categoricalvariablesidx))
		{ 
			cat("Note: categorical variables are found in data. Please use option categoricalvariablesidx = 'all' to match them more closely.\n")
		}		 
		if (is.data.frame(X)) 
		{ 	cat("X is a data frame. String or factors have been converted to numeric values.\n") }
		X <- NAfactor2matrix(X, toGrep = "anythingParticular")
		if (!is.null(categoricalvariablesidx))
		{
			if (categoricalvariablesidx[1] == "all")
			{	
				factorVariables <- which(getFactors > 0)
				if (length(factorVariables) > 0) 	{ 	categoricalvariablesidx = factorVariables	}
				else {  cat("\nNo categorical variables found. Please type them manually\n")  }
			}
			else 
			{  
				if (is.character(categoricalvariablesidx[1]))
				{	categoricalvariablesidx = sort(match(categoricalvariablesidx, colnames(X)))	}
			}
		}
		if ( length(which( (X == Inf) | (X == -Inf) ) ) > 0) 
		{ stop("Inf or -Inf values found in data. Learning can not be done.\nRemove or replace them with NA in order to learn.\n") }
		XY <- NATreatment(X, Y, na.action = na.action, regression = regression)
		X  <- XY$X
		Y  <- XY$Y
		rm(XY)		
		if (randomcombination[1] > 0)
		{ 
			randomcombinationString = randomcombination[1]
			L.combination = length(randomcombination)
			if (L.combination%%3 != 0)
			{	weights = round(replicate(length(randomcombination)/2, sample(c(-1,1),1)*runif(1)), 2)	}
			else
			{ 	
				weights = round(randomcombination[(L.combination + 1 - L.combination/3):L.combination],2)
				randomcombination = randomcombination[1:(L.combination - (L.combination/3))]
			}			
			for (i in 2:length(randomcombination)) {  randomcombinationString = paste(randomcombinationString, randomcombination[i], sep=",") }
			for (i in 1:length(weights)) {  randomcombinationString = paste(randomcombinationString, weights[i], sep=",") }
			if (!is.null(xtest)) { xtest <- randomCombination(NAfactor2matrix(xtest, toGrep = "anythingParticular"), combination = randomcombination, weights = weights) }			
			randomcombinationObject = list(randomcombination, weights)					
			X <- randomCombination(X, combination = randomcombinationObject[[1]], weights = randomcombinationObject[[2]])			
		}
		else
		{   randomcombinationObject = randomcombination }
	}	
	RUF.model <- randomUniformForestCore(X, trainLabels = Y, ntree = ntree, nodeMinSize = nodesize, maxNodes = maxnodes, 
		features = mtry, rf.bootstrap = replace, depth = depth, depthControl = depthcontrol, 
		rf.treeSubsampleRate = subsamplerate, classwt = classwt, classCutOff = classcutoff, 
		rf.overSampling = oversampling, rf.targetClass = targetclass, rf.rebalancedSampling = rebalancedsampling, 
		rf.outputPerturbationSampling = outputperturbationsampling,	rf.randomCombination = randomcombinationObject, 
		rf.randomFeature = randomfeature, rf.treeBagging = bagging, rf.featureSelectionRule = featureselectionrule[1], 
		rf.regression = regression, use.OOB = OOB, BreimanBounds = BreimanBounds, variableImportance = importance, 
		whichCatVariables = categoricalvariablesidx, logX = logX, threads = threads, useSubTrees = usesubtrees,unsupervised = unsupervised,
		parallelPackage = parallelpackage[1])
	if (!is.null(classwt))
	{ 
		classwtString = classwt[1]
		for (i in 2:length(classwt)) { 	classwtString = paste(classwtString, classwt[i], sep=",") }	
	}	
	if (length(targetclass) > 1)
	{ 
		targetclassString = targetclass[1]
		for (i in 2:length(targetclass)) { 	targetclassString = paste(targetclassString, targetclass[i], sep=",") }	
	}	
	if (length(rebalancedsampling) > 1)
	{ 
		rebalancedsamplingString = rebalancedsampling[1]
		for (i in 2:length(rebalancedsampling)) 
		{ rebalancedsamplingString = paste(rebalancedsamplingString, rebalancedsampling[i], sep= ",") }	
	}	
	if (RUF.model$regression) {	classcutoff = c(0,0)  }	
	if (as.numeric(classcutoff[2]) == 0) {	classcutoffString = FALSE }
	else
	{
		classcutoffString = levels(Y)[which(levels(Y) == as.character(classcutoff[1]))]
		classcutoffString = paste("Class ", classcutoffString, "," , as.numeric(classcutoff[2])*100, "%", sep ="")
	}		
	paramsObject = c(ntree, mtry, nodesize, maxnodes, as.character(replace), bagging, depth, 
		ifelse(is.null(depthcontrol), length(depthcontrol)> 0, depthcontrol), OOB, importance, subsamplerate, 
		ifelse(is.null(classwt), length(classwt) > 0, classwtString), classcutoffString, 
		ifelse((oversampling == 0), (oversampling != 0), oversampling), outputperturbationsampling, 
		ifelse(length(targetclass) > 1, targetclassString, targetclass), 
		ifelse(length(rebalancedsampling) > 1, rebalancedsamplingString, rebalancedsampling), 
		ifelse((randomcombination[1] == 0),(randomcombination != 0), randomcombinationString), randomfeature, 
		ifelse(is.null(categoricalvariablesidx), length(categoricalvariablesidx) > 0, length(categoricalvariablesidx)), featureselectionrule[1])				
	if (RUF.model$regression) 
	{ 
		paramsObject[21] = if (featureselectionrule[1] == "L1") { "Sum of absolute residuals" }  else { "Sum of squared residuals" }
	    if ((paramsObject[5] == "FALSE") & (subsamplerate == 1)) { paramsObject[8] = FALSE }
	}	
	names(paramsObject) = c("ntree", "mtry", "nodesize", "maxnodes", "replace", "bagging", "depth", "depthcontrol", "OOB", "importance", "subsamplerate", "classwt", "classcutoff", "oversampling", "outputperturbationsampling", "targetclass", "rebalancedsampling", "randomcombination", "randomfeature", "categorical variables", "featureselectionrule")						
	paramsObject = as.data.frame(paramsObject)	
	RUF.model$logX = logX	
	if (!regression & !is.null(ytest) ) 
	{ 
		if (!is.factor(ytest))
		{ ytest = as.factor(ytest) }
	}
	RUFObject <- genericOutput(xtest, ytest, paramsObject, RUF.model, ytrain = Y, classcutoff = classcutoff)
	RUFObject$logX = logX
	if (!is.null(Y)) { 	RUFObject$y = Y	 }	
	if (is.null(colnames(X))) 
	{  
		varNames = NULL
		for (i in 1:ncol(X)) { varNames = c(varNames, paste("V", i, sep="")) }
		RUFObject$variablesNames = varNames
	}
	else
	{	RUFObject$variablesNames = colnames(X)	}	
	if (randomcombination[1] > 0)
	{ 	
		tempVarNames = vector(length = length(randomcombination)/2)
		idx = 1
		for (i in 1:(length(randomcombination)/2)) 
		{ 
			tempVarNames[i] = paste("V", randomcombination[idx], "x", randomcombination[idx+1], sep="") 
			idx = idx + 2
		}
		RUFObject$variablesNames = c(RUFObject$variablesNames, tempVarNames)
	}
	if (!is.null(categoricalvariablesidx)) {  RUFObject$categoricalvariables = categoricalvariablesidx }
	if (unsupervised) 
	{ 
		RUFObject$unsupervised = TRUE 
		RUFObject$unsupervisedMethod  = if (length(unsupervisedMethod) == 3) {unsupervisedMethod[1] } 
							else { if (length(unsupervisedMethod) == 1) { unsupervisedMethod } else { unsupervisedMethod[1:2] } }
	}
	RUFObject$call <- match.call()
	class(RUFObject) <- "randomUniformForest"
	RUFObject
}

randomUniformForestCore <- function(trainData, trainLabels = 0, 
	features = floor(4/3*(ncol(trainData))), 
	ntree = 100, 
	nodeMinSize = 1,
	maxNodes = Inf,
	depth = Inf,
    depthControl = NULL,
	splitAt = "random", 
	rf.regression = TRUE, 
	rf.bootstrap = ifelse(rf.regression, FALSE, TRUE), 
	use.OOB = TRUE,
	BreimanBounds = ifelse(use.OOB, TRUE, FALSE),
	rf.treeSubsampleRate = ifelse(rf.regression, 0.7, 1),
	rf.treeBagging = FALSE, 
	classwt = NULL,
	classCutOff = c(0,0),
	rf.overSampling = 0,
	rf.targetClass = -1,
	rf.outputPerturbationSampling = FALSE,
    rf.rebalancedSampling = FALSE,
	rf.featureSelectionRule = c("entropy", "gini", "random", "L2", "L1"),
	variableImportance = TRUE,	
	rf.randomCombination = 0,
	rf.randomFeature = FALSE,
	whichCatVariables = NULL,
	logX = FALSE,
	unsupervised = FALSE,
	useSubTrees = FALSE,
	threads = "auto",	
	parallelPackage = "doParallel",
	export = c("uniformDecisionTree", "CheckSameValuesInAllAttributes", "CheckSameValuesInLabels", "fullNode", "genericNode", "leafNode", "randomUniformForestCore.predict", "onlineClassify", "overSampling", "predictDecisionTree", "options.filter", "majorityClass", "randomCombination", "randomWhichMax", "which.is.na", "factor2vector", "outputPerturbationSampling", "rmNA", "count.factor", "find.idx", "classifyMatrixCPP", 
	"L2DistCPP", "checkUniqueObsCPP", "crossEntropyCPP", "giniCPP", "L2InformationGainCPP", 
	"entropyInformationGainCPP", "runifMatrixCPP", "NATreatment", "rmInf", "rm.InAList"), 
	...)
{
	set.seed(sample(ntree,1))
	n = nrow(trainData)
	p = ncol(trainData)
	if (exists("categorical", inherits = FALSE)) { whichCatVariables = categorical }
	else { categorical = NULL  }
	if (useSubTrees) 
	{ 
		if (depth == Inf) { depth = min(10, floor(0.5*log(n)/log(2))) }
		if (rf.outputPerturbationSampling)
		{
			rf.outputPerturbationSampling = FALSE
			cat("'output perturbation sampling' is not currently compatible with 'usesubtrees'.\n Option has been reset.\n")
		}
	}
	if ((rf.randomFeature) & (features == "random") ) { stop("random feature is a special case of mtry = 'random'") }
	if (is.numeric(threads))
	{
		if (threads < 1) { stop("Number of threads must be positive") }
	}
	if (!is.matrix(trainData)) { trainData <- NAfactor2matrix(trainData, toGrep = "anythingParticular") }
	if (!is.null(whichCatVariables)) 
	{ 
	  cat("Considering use of Dummies (for example by using formula) for categorical variables might be an alternative.\n")
	  # rmFromCatIdx = apply(trainData[,whichCatVariables], 2,function(Z) length(table(Z)))
	  # whichCatVariables = whichCatVariables[which(rmFromCatIdx != 2)]
	  # if (length(whichCatVariables) == 0) {  whichCatVariables = NULL }
	}
	{
		if (!is.null(classwt)) 
		{ 
			classwt = classwt/sum(classwt) 
			if (is.na(sum(classwt))) { stop ("NA found in class weights.\n") }			
			if (sum(classwt) > 1.1)  { stop ("(inverse) weights do not sum to 1.\n") }
		}		
		if (length(trainLabels) != 1)
		{	
			if (n != length(trainLabels)) 
			{ stop ("X and Y don't have the same size.\n")	}	
		}		
		if (!is.matrix(trainData))
		{ stop ("X cannot be converted to a matrix. Please provide true matrix not data frame.") }		
		if (rf.treeSubsampleRate == 0) { rf.treeSubsampleRate = 1 }		
		if (rf.treeBagging) { features = min(features, p)  }
		if ( (rf.treeSubsampleRate < 0.5) & (n < 300) ) 
		{ rf.treeSubsampleRate = 0.5;  cat("Too small output matrix. Subsample rate has been set to 0.5.\n") }		
		if (!is.character(nodeMinSize))
		{
			if ( (nodeMinSize != 1) & (nodeMinSize >  floor(nrow(trainData)/4)) ) { stop("nodeMinSize is too high. Not suitable for a random forest.\n") }	
			if ( (nodeMinSize < 1) ) {  nodeMinSize  = 1; cat("Minimal node size has been set to 1.\n") }
		}		
		if ( features == 1 ) { rf.randomFeature = TRUE }
		if ( features < 1 )  
		{ 
			features = floor(4/3*p)  
			bagging = FALSE
			cat("Error setting mtry. Resetting to default values.\n") 
		}			
		if (is.factor(trainLabels))
		{ 	
			if ((as.numeric(classCutOff[2]) != 0)) 
			{  	
				classCutOff = c(which(levels(trainLabels) == as.character(classCutOff[1])), 
				0.5/as.numeric(classCutOff[2])) 
				if (length(classCutOff) == 1)
				{ stop("Label not found. Please provide name of the label instead of its position.\n") }
			}					
			rf.regression = FALSE 
			labelsObject = factor2vector(trainLabels)
			trainLabels = as.numeric(as.vector(labelsObject$vector))
			if (length(unique(trainLabels)) == 1 ) { stop ("Y is a constant value. Thus, learning is not needed.\n")	}			
		}
		else
		{	
			if (rf.regression & (length(unique(trainLabels)) > 32)) 
			{ 
				if ((rf.treeSubsampleRate < 1) & (rf.bootstrap == TRUE))
				{ cat("For only accuracy, use option 'subsamplerate = 1' and 'replace = FALSE'\n") }				
				classCutOff =  c(0,0) 				
			}
			else
			{
				if (!rf.regression)
				{
					if ((as.numeric(classCutOff[2]) != 0)) 
					{  	
						classCutOff = c(which(levels(trainLabels) == as.character(classCutOff[1])), 0.5/as.numeric(classCutOff[2])) 
						if (length(classCutOff) == 1)
						{ stop("Label not found. Please provide name of the label instead of its position") }
					}					
					labelsObject = factor2vector(trainLabels)
					trainLabels = as.numeric(as.vector(labelsObject$vector))
					labelsObject$factors = levels(as.factor(trainLabels))
					if (length(unique(trainLabels)) == 1 ) { stop ("Y is a constant value. Thus, learning is not needed\n")	}
				}
				else
				{		
					if (length(unique(trainLabels)) <= 32)  
					{ 
						cat("Regression has been performed but there is less than 32 distinct values.\nRegression option could be set to FALSE, if classification is #required.\nIf outputs are factors, conversion is automatically done.\n")
						classCutOff = c(0,0)
					}
				}
			}
		}							
		if (!rf.regression)
		{
			rf.classes <- as.numeric(levels(as.factor(trainLabels)) )			
			if (as.character(rf.classes[1]) != labelsObject$factors[1]) 
			{	
				cat("Labels", labelsObject$factors, "have been converted to", rf.classes, 
				"for ease of computation and will be used internally as a replacement.\n")	
			}			
			if (!is.null(classwt))
			{ 	
				if (length(classwt) != length(rf.classes)) 
				{ stop("Length of class weights is not equal to length of classes.\n") }
			}
		}
		else
		{ 	rf.classes = trainLabels; rf.overSampling = 0; rf.rebalancedSampling = FALSE; classCutOff = c(0,0)  }		
		if (logX)
		{ 
		   if (is.null(whichCatVariables))  {  trainData <- generic.log(trainData) }
		   else   {  trainData[,-whichCatVariables] <- generic.log(trainData[,-whichCatVariables]) }
		}		
		if (!rf.regression)  
		{ 
			if ( (rf.featureSelectionRule[1] == "L1") | (rf.featureSelectionRule[1] == "L2") )
			{	
				rf.featureSelectionRule[1] = "entropy" 
				cat("Feature selection rule has been set to entropy.\n")
			}
		}
		else 
		{  
			if ( (rf.featureSelectionRule[1] == "entropy")  |  (rf.featureSelectionRule[1] == "gini") ) 
			{	rf.featureSelectionRule[1] = "L2" }			
			if (!is.null(depthControl) & (!is.character(depthControl)))
			{	
				if (depthControl < 1)
				{
					depthControl = NULL
					cat("'depthcontrol' option is lower than its lower bound. Resetting to default value.\n")
				}				
			}			
		}
		if (!rf.bootstrap & (rf.treeSubsampleRate == 1)) {  use.OOB = FALSE }
	}
	{
		## #require(parallel)	
		max_threads = detectCores()		
		if (threads == "auto")
		{	
			if (max_threads == 2) { threads = max_threads }
			else {	threads  = max(1, max_threads - 1)  }
		}
		else
		{
			if (max_threads < threads) 
			{	cat("Note: number of threads is higher than logical threads in this computer.\n") }
		}
		
		{
			## #require(doParallel)			
			Cl = makePSOCKcluster(threads, type = "SOCK")
			registerDoParallel(Cl)
		}
		chunkSize  <-  ceiling(ntree/getDoParWorkers())
		smpopts  <- list(chunkSize = chunkSize)
	}
	rufObject = vector('list', ntree)
	if (use.OOB)
	{
		rufObject <- foreach(i = 1:ntree, .export = export, .options.smp = smpopts, .inorder = FALSE, 
		.multicombine = TRUE, .maxcombine = ntree) %dopar%
		{
			uniformDecisionTree(trainData, trainLabels, nodeMinSize = nodeMinSize, maxNodes = maxNodes, 
			treeFeatures = features, getSplitAt = splitAt, regression = rf.regression, bootstrap = rf.bootstrap, treeSubsampleRate = rf.treeSubsampleRate, treeDepth = depth, treeClasswt = classwt,
			treeOverSampling = rf.overSampling, targetClass = rf.targetClass, OOB = use.OOB, 
			treeRebalancedSampling = rf.rebalancedSampling, treeBagging = rf.treeBagging, 
			randomFeature = rf.randomFeature, treeCatVariables = whichCatVariables, 
			outputPerturbation = rf.outputPerturbationSampling, featureSelectionRule = rf.featureSelectionRule, treeDepthControl = depthControl, unsupervised = unsupervised)
		}		
		stopCluster(Cl)
		if (useSubTrees)
		{
			rufObject2 = vector('list', ntree)
			Cl = makePSOCKcluster(threads, type = "SOCK") 
			registerDoParallel(Cl)
						
			rufObject2 <- foreach(i = 1:ntree, .export = export, .options.smp = smpopts, .inorder = FALSE, 
			.multicombine = TRUE, .maxcombine = ntree) %dopar%
			{
				finalTree = list()
				OOBIdx = rufObject[[i]]$OOB.idx
				idxList = rufObject[[i]]$idxList
				trainData2 = trainData[rufObject[[i]]$followIdx,]
				trainLabels2 = trainLabels[rufObject[[i]]$followIdx]
				nodeVector = rufObject[[i]]$nodeVector
				lengthIdx = sapply(idxList, length)
				rmlengthIdx = which(lengthIdx < 2)
				nodeVector = if (length(rmlengthIdx) > 0) nodeVector[-rmlengthIdx] else nodeVector
				newIdxList <- if (length(rmlengthIdx) > 0) rm.InAList(idxList, rmlengthIdx) else idxList
				originalTree = rufObject[[i]]$Tree
				lengthNodeVector = length(nodeVector)
				newTrees = vector('list', lengthNodeVector)
								
				if ( (lengthNodeVector > 0) )
				{
					newTrees <- lapply(newIdxList, function(Z)
						{
							newTrainData = trainData2[Z, ,drop = FALSE]; newTrainLabels = trainLabels2[Z]
													
							uniformDecisionTree(newTrainData, newTrainLabels, nodeMinSize = nodeMinSize, 
							maxNodes = maxNodes, treeFeatures = features, getSplitAt = splitAt, 
							regression = rf.regression, bootstrap = FALSE, treeSubsampleRate = 1, treeDepth = Inf, treeClasswt = classwt, treeOverSampling = rf.overSampling, targetClass = rf.targetClass, OOB = FALSE, treeRebalancedSampling = rf.rebalancedSampling, treeBagging = rf.treeBagging, 
							randomFeature = rf.randomFeature, treeCatVariables = whichCatVariables, 
							outputPerturbation = rf.outputPerturbationSampling, 
							featureSelectionRule = rf.featureSelectionRule, treeDepthControl = NULL, 
							unsupervised = unsupervised, moreThreadsProcessing = FALSE, moreNodes = TRUE)
						}
					)
					
					for (j in 1:lengthNodeVector)
					{
						toInsertTree = newTrees[[j]]$Tree 
						keepMaxIdx = nrow(originalTree)
						toInsertTree[,1:2] = toInsertTree[,1:2] + keepMaxIdx
						toInsertTree[toInsertTree[,"status"] == -1,1:2] = 0
						originalTree = rbind(originalTree, toInsertTree)
						originalTree[nodeVector[j],] = toInsertTree[1,]
					}
					rownames(originalTree) = 1:nrow(originalTree)
				}
				finalTree$Tree = originalTree
				finalTree$OOB.idx = OOBIdx
				finalTree
			}
			rufObject = rufObject2 
			stopCluster(Cl)
		}
		{
			OOB.matrix = NULL 
			new.rufObject = vector("list", ntree)
			new.rufObject$Tree = lapply(rufObject, function(Z) Z$Tree) 
			for (i in 1:ntree) 	
			{	OOB.matrix <- rbind(OOB.matrix, cbind(rufObject[[i]]$OOB.idx, rep(i,length(rufObject[[i]]$OOB.idx))))  } 
			OOB.val = sort(unique(OOB.matrix[,1])) 
			n.OOB = nrow(trainData) 
			OOB.votes2 = matrix(Inf, n.OOB, ntree)			
			if (is.null(classwt))
			{
				OOB.votes <- randomUniformForestCore.predict(new.rufObject$Tree, trainData, 
					pr.regression = rf.regression, classes = rf.classes, OOB.idx = TRUE, 
					pr.parallelPackage = parallelPackage[1], pr.imbalance = classCutOff, pr.threads = threads)			 
				for (j in 1:ntree)
				{
					idxJ = which(OOB.matrix[,2] == j)
					if (length(idxJ) > 0) {	OOB.votes2[OOB.matrix[idxJ,1],j] = OOB.votes[OOB.matrix[idxJ,1],j] }
				}
				OOB.object <- majorityClass(OOB.votes2, rf.classes, m.imbalance = classCutOff, 
					m.regression = rf.regression)
			}
			else
			{
				OOB.allWeightedVotes = matrix(Inf, n.OOB, ntree)
				OOB.votes <- randomUniformForestCore.predict(new.rufObject$Tree, trainData, 
					pr.regression = rf.regression, classes = rf.classes, pr.classwt = TRUE, OOB.idx = TRUE, pr.parallelPackage = parallelPackage[1], pr.imbalance = classCutOff, pr.threads = threads)			
				for (j in 1:ntree)
				{
					idxJ = which(OOB.matrix[,2] == j)
					if (length(idxJ) > 0) 	
					{	
						OOB.votes2[OOB.matrix[idxJ,1],j] = OOB.votes$all.votes[OOB.matrix[idxJ,1],j] 	
						OOB.allWeightedVotes[OOB.matrix[idxJ,1],j] = OOB.votes$allWeightedVotes[OOB.matrix[idxJ,1],j]
					}
				}				
				OOB.object <- majorityClass(OOB.votes2, rf.classes, m.classwt = OOB.allWeightedVotes, 
					m.imbalance = classCutOff, m.regression = rf.regression)			
			}			
			OOB.pred = OOB.object$majority.vote
			if (BreimanBounds)
			{
				strengthCorr.object <- strength_and_correlation(s.trainLabels = trainLabels, OOB.votes2, 
					OOB.object, rf.classes, s.regression = rf.regression, s.parallelPackage = parallelPackage[1], s.threads = threads)
			}
			if (rf.regression)
			{  
				OOB.confusion = "only prediction error for regression. See ..$pred.error"				
				MSE <- L2Dist(OOB.pred[OOB.val], trainLabels[OOB.val])/n.OOB				
				pred.error = round(MSE,4)
				percent.varExplained = max(0, 100*round(1 - MSE/var(trainLabels[OOB.val]),4))
			}
			else
			{
				if ( min(trainLabels) == 0) { OOB.pred = OOB.pred - 1 }			
				if ( (length(unique(trainLabels[OOB.val])) == 1) | (length(unique(OOB.pred[OOB.val])) == 1) )
				{	
					OOB.confusion = "Minority class can not be predicted or OOB predictions contain only one label."
					pred.error = "Minority class can not be predicted or OOB predictions contain only one label."
				}
				else
				{
					OOB.confusion <- confusion.matrix(OOB.pred[OOB.val], trainLabels[OOB.val])
					pred.error <- generalization.error(OOB.confusion)
				}
			}
		}
		if (variableImportance)
		{	
			Cl = makePSOCKcluster(threads, type = "SOCK")
			registerDoParallel(Cl)
			gen.rufObject = new.rufObject$Tree	
		}
		else
		{	
			if (rf.regression)
			{ 
				if (BreimanBounds)
				{
					return(list(object = new.rufObject$Tree, OOB = OOB.confusion, OOB.predicts = as.numeric(OOB.pred), 
					OOB.strengthCorr = strengthCorr.object, OOB.votes = OOB.votes2,	pred.error = pred.error, 
					percent.varExplained = percent.varExplained, regression = rf.regression) )
				}
				else
				{
					return(list(object = new.rufObject$Tree, OOB = OOB.confusion, OOB.predicts = as.numeric(OOB.pred), 
					OOB.votes = OOB.votes2,	pred.error = pred.error, percent.varExplained = percent.varExplained,  
					regression = rf.regression) )
				}
			}
			else
			{ 	
				if (BreimanBounds)
				{
					return(list(object = new.rufObject$Tree, OOB = OOB.confusion, OOB.predicts = as.numeric(OOB.pred), 
						OOB.strengthCorr = strengthCorr.object, OOB.votes = OOB.votes2, pred.error = pred.error, 
						regression = rf.regression) )
				}
				else
				{
					return(list(object = new.rufObject$Tree, OOB = OOB.confusion, OOB.predicts = as.numeric(OOB.pred), 
						OOB.votes = OOB.votes2, pred.error = pred.error, regression = rf.regression) )
				}					
			}
		}		
	} 
	else
	{
		if (!useSubTrees)
		{
			rufObject <- foreach(i = 1:ntree, .export = export, .options.smp  = smpopts, .inorder = FALSE, 
			.multicombine = TRUE, .maxcombine = ntree) %dopar%	
			{	
				uniformDecisionTree(trainData, trainLabels, nodeMinSize = nodeMinSize, maxNodes = maxNodes, 
				treeFeatures = features, getSplitAt = splitAt, regression = rf.regression, bootstrap = rf.bootstrap, treeSubsampleRate = rf.treeSubsampleRate, treeDepth = depth, treeClasswt = classwt, 
				treeOverSampling = rf.overSampling, targetClass = rf.targetClass, treeBagging = rf.treeBagging, 
				treeRebalancedSampling = rf.rebalancedSampling, randomFeature = rf.randomFeature, 
				treeCatVariables = whichCatVariables, outputPerturbation = rf.outputPerturbationSampling, 
				featureSelectionRule = rf.featureSelectionRule, treeDepthControl = depthControl, 
				unsupervised = unsupervised, moreThreadsProcessing = useSubTrees)$Tree
			}
		}		
		else
		{
			rufObject <- foreach(i = 1:ntree, .export = export, .options.smp  = smpopts, .inorder = FALSE, 
			.multicombine = TRUE, .maxcombine = ntree) %dopar%	
			{	
				uniformDecisionTree(trainData, trainLabels, nodeMinSize = nodeMinSize, maxNodes = maxNodes, 
				treeFeatures = features, getSplitAt = splitAt, regression = rf.regression, bootstrap = rf.bootstrap, treeSubsampleRate = rf.treeSubsampleRate, treeDepth = depth, treeClasswt = classwt, 
				treeOverSampling = rf.overSampling, targetClass = rf.targetClass, treeBagging = rf.treeBagging, 
				treeRebalancedSampling = rf.rebalancedSampling, randomFeature = rf.randomFeature, 
				treeCatVariables = whichCatVariables, outputPerturbation = rf.outputPerturbationSampling, 
				featureSelectionRule = rf.featureSelectionRule, treeDepthControl = depthControl, 
				unsupervised = unsupervised, moreThreadsProcessing = useSubTrees)
			}
						
			stopCluster(Cl)
			rufObject2 = vector('list', ntree)
			Cl = makePSOCKcluster(threads,  type = "SOCK") 
			registerDoParallel(Cl)
						
			rufObject2 <- foreach(i = 1:ntree, .export = export, .options.smp = smpopts, .inorder = FALSE, 
			.multicombine = TRUE, .maxcombine = ntree) %dopar%
			{
				OOBIdx = rufObject[[i]]$OOB.idx
				idxList = rufObject[[i]]$idxList
				trainData2 = trainData[rufObject[[i]]$followIdx,]
				trainLabels2 = trainLabels[rufObject[[i]]$followIdx]
				nodeVector = rufObject[[i]]$nodeVector
				lengthIdx = sapply(idxList, length)
				rmlengthIdx = which(lengthIdx < 2)
				nodeVector = if (length(rmlengthIdx) > 0) nodeVector[-rmlengthIdx] else nodeVector
				newIdxList <- if (length(rmlengthIdx) > 0) rm.InAList(idxList, rmlengthIdx) else idxList
				originalTree = rufObject[[i]]$Tree
				lengthNodeVector = length(nodeVector)
				newTrees = vector('list', lengthNodeVector)
				
				if (length(nodeVector) > 0)
				{
					newTrees <- lapply(newIdxList, function(Z)
						{
							newTrainData = trainData2[Z, ,drop = FALSE]; newTrainLabels = trainLabels2[Z]
							
							uniformDecisionTree(newTrainData, newTrainLabels, nodeMinSize = nodeMinSize, 
							maxNodes = maxNodes, treeFeatures = features, getSplitAt = splitAt, 
							regression = rf.regression, bootstrap = FALSE, treeSubsampleRate = 1, treeDepth = Inf, treeClasswt = classwt, treeOverSampling = rf.overSampling, targetClass = rf.targetClass, OOB = FALSE, treeRebalancedSampling = rf.rebalancedSampling, treeBagging = rf.treeBagging, 
							randomFeature = rf.randomFeature, treeCatVariables = whichCatVariables, 
							outputPerturbation = rf.outputPerturbationSampling,
							featureSelectionRule = rf.featureSelectionRule, treeDepthControl = NULL, 
							unsupervised = unsupervised, moreThreadsProcessing = FALSE, moreNodes = TRUE)
						}
					)
					
					for (j in 1:lengthNodeVector)
					{
						toInsertTree = newTrees[[j]]$Tree 
						keepMaxIdx = nrow(originalTree)
						toInsertTree[,1:2] = toInsertTree[,1:2] + keepMaxIdx
						toInsertTree[toInsertTree[,"status"] == -1,1:2] = 0
						originalTree = rbind(originalTree, toInsertTree)
						originalTree[nodeVector[j],] = toInsertTree[1,]
					}
					rownames(originalTree) = 1:nrow(originalTree)
				}				
				originalTree
			}
			rufObject = rufObject2 
		}
					
		if (variableImportance) {  gen.rufObject = rufObject	}
		else 
		{	
			stopCluster(Cl) 			
			return( list(object = rufObject, regression = rf.regression))	
		}
	}
	if (variableImportance)
	{
	   	if (ntree <= 100) {  threads = 1 }						
		varImpMatrix1 <- unlist(lapply(gen.rufObject, function(Z) Z[,"split var"]))		
		if (rf.regression) 	{  varImpMatrix2 <- unlist(lapply(gen.rufObject, function(Z) Z[,"L2Dist"]))/1000 }
		else {	varImpMatrix2 <- unlist(lapply(gen.rufObject, function(Z) Z[,"Gain"])) }			
		varImpMatrix <- cbind(varImpMatrix1, varImpMatrix2)
		if (!rf.regression)
		{   
			predMatrix <- foreach(i = 1:ntree, .options.smp = smpopts, .inorder = TRUE, .combine = rbind, 
			.multicombine = TRUE) %dopar%	
			{ 	
				predIdx = which(gen.rufObject[[i]][,"left daughter"] == 0)
				predClass = gen.rufObject[[i]][predIdx, "prediction"]				
				predVar = vector(length = length(predIdx))
				if (useSubTrees)
				{
					for (j in seq_along(predIdx))
					{ 
						predVar[j] = gen.rufObject[[i]][which(gen.rufObject[[i]][,"left daughter"] == predIdx[j] |gen.rufObject[[i]][,"right daughter"] == predIdx[j]), "split var"][1]
					}
				}
				else
				{
					for (j in seq_along(predIdx))
					{ 
						if ((predIdx[j] %% 2) == 0) 
						{ 
							predVar[j] = gen.rufObject[[i]][which(gen.rufObject[[i]][,"left daughter"] == predIdx[j]), "split var"] 
						}
						else 
						{	
							predVar[j] = gen.rufObject[[i]][which(gen.rufObject[[i]][,"right daughter"] == predIdx[j]), "split var"] 
						}
					}
				}			
				cbind(predVar, predClass)
			}	
		}
		stopCluster(Cl)
		varImpMatrix = varImpMatrix[-which(varImpMatrix[,1] == 0),]		
		na.gain = which(is.na(varImpMatrix[,2]))
		if (length(na.gain) > 0) { 	varImpMatrix = varImpMatrix[-na.gain,]	}			
		rf.var = unique(sortCPP(varImpMatrix[,1]))
		n.var = length(rf.var)			
		if (!rf.regression)
		{	
			var.imp.output <- matrix(NA, n.var, 4)
			for (i in 1:n.var)
			{  	
				classTable = sort(table(predMatrix[which(rf.var[i] == predMatrix[,1]),2]),decreasing = TRUE)
				classMax = as.numeric(names(classTable))[1]
				classFreq =	(classTable/sum(classTable))[1]
				var.imp.output[i,] = c(rf.var[i], sum(varImpMatrix[which(rf.var[i] == varImpMatrix[,1]),2]), 
					classMax, classFreq)
			}
		}
		else
		{	
			var.imp.output <- matrix(NA, n.var, 2)
			for (i in 1:n.var)
			{  	var.imp.output[i,] = c(rf.var[i], sum(varImpMatrix[which(rf.var[i] == varImpMatrix[,1]),2]))	}
		}
		var.imp.output = var.imp.output[order(var.imp.output[,2], decreasing = TRUE),]
		var.imp.output = round(cbind( var.imp.output, 100*var.imp.output[,2]/max(var.imp.output[,2])),2)
		percent.importance = round(100*var.imp.output[,2]/sum(var.imp.output[,2]),0)
		var.imp.output[,2] = round(var.imp.output[,2],0)		
		var.imp.output = cbind(var.imp.output, percent.importance)		
		if (rf.regression) 
		{	colnames(var.imp.output) = c("variables", "score", "percent", "percent importance")	}
		else
		{	
			colnames(var.imp.output) = c("variables", "score", "class", "class frequency", "percent", 
				"percent importance")
			row.names(var.imp.output) = NULL	
		}		
		var.imp.output = data.frame(var.imp.output)		
		if (!is.null(colnames(trainData))) 	{	var.imp.output[,1] = colnames(trainData)[var.imp.output[,1]]	}		
		if (use.OOB)
		{	
			if (rf.regression)
			{
				if (BreimanBounds)
				{
					return(list(object = gen.rufObject, OOB = OOB.confusion, OOB.predicts = as.numeric(OOB.pred), 
					OOB.strengthCorr = strengthCorr.object, OOB.votes = OOB.votes2, pred.error = pred.error, 
					percent.varExplained = percent.varExplained,  variableImportance = var.imp.output, 
					regression = rf.regression) )
				}
				else
				{
					return(list(object = gen.rufObject, OOB = OOB.confusion, OOB.predicts = as.numeric(OOB.pred), 
					OOB.votes = OOB.votes2, pred.error = pred.error, percent.varExplained = percent.varExplained, 
					variableImportance = var.imp.output, regression = rf.regression) )
				}
			}
			else
			{
				if (BreimanBounds)
				{
					return( list(object = gen.rufObject, OOB = OOB.confusion, OOB.predicts = as.numeric(OOB.pred), 
					OOB.strengthCorr = strengthCorr.object, OOB.votes = OOB.votes2, pred.error = pred.error, 
					variableImportance = var.imp.output, regression = rf.regression) )
				}
				else
				{
					return(list(object = gen.rufObject, OOB = OOB.confusion, OOB.predicts = as.numeric(OOB.pred), 
					OOB.votes = OOB.votes2, pred.error = pred.error, variableImportance = var.imp.output, regression = rf.regression) )	
				}
			}
		}
		else
		{	return(list(object = gen.rufObject, variableImportance = var.imp.output, regression = rf.regression) )	}
	}
}

rUniformForestPredict <- function(object, X, 
	type = c("response", "prob", "votes", "confInt", "ranking", "quantile", "truemajority", "all"),
	classcutoff = c(0,0), 
	conf = 0.95,
	whichQuantile = NULL,
	rankingIDs = NULL,
	threads = "auto", 
	parallelpackage = "doParallel", ...)
{
	object <- filter.object(object)
	X <- if (is.vector(X)) { t(X) } else  { X }
	X <- fillVariablesNames(X)	
	if (!is.null(object$formula) & (length(object$variablesNames) != ncol(X)))
	{
		mf <- model.frame(data = as.data.frame(cbind(rep(1, nrow(X)),X)))
		X <- model.matrix(attr(mf, "terms"), data = mf)[,-1]
	}
	else
	{
		if (length(object$variablesNames) != ncol(X)) 
		{	
			cat("Data to predict have not the same dimension than data that have been computed by the model.\n Relevant variables have been extracted.\n")
			newVars <- rmNA(match(object$variablesNames, colnames(X)))
			if (length(newVars) == 0) 
			{ 
				stop("No relevant variable has been found. Please give to both train and test data the same column names.\n") 
			}
			X = X[,newVars, drop = FALSE]
		}
	}
	classes = object$classes	
	flag1 = flag2 = 1
	for (i in 1:ncol(object$forestParams))
	{			
		classCutOffString = as.vector(object$forestParams[which(row.names(object$forestParams) == "classcutoff"), i])		
		if ((classCutOffString != "FALSE") & (as.numeric(classcutoff[2]) == 0))
		{
			classCutOffString = rm.string(classCutOffString, "%")
			classCutOffString = rm.string(classCutOffString, "Class ")
			classCutOffString = strsplit(classCutOffString, ",")
			classcutoff[1] = which(classes == classCutOffString[[1]][1])
			classcutoff[2] = 0.5/(as.numeric(classCutOffString[[1]][2])/100)
		}
		else
		{	
			if ((as.numeric(classcutoff[2]) != 0) & (i == 1)) 
			{  	
				classcutoff <- c(which(classes == as.character(classcutoff[1])), 0.5/(as.numeric(classcutoff[2]))) 
				if (length(classcutoff) == 1) 
				{ stop("Label not found. Please provide name of the label instead of its position") }
				if (i > 1) { cat("Only last cutoff will be used in incremental random Uniform Forest.\n") }
			}		
		}				
		classwtString = as.vector(object$forestParams[which(row.names(object$forestParams) == "classwt"),i])	
		classwt = FALSE
		if (is.na(as.logical(classwtString))) {  classwt = TRUE; flag1 = flag2 = -1  }
		else 
		{	flag2 = 1	}		
		if (flag1 != flag2) 
		{ stop("Incremental random Uniform forest. Class reweighing must remain (or miss) in all forests.\n") }	
		randomCombinationString = as.vector(object$forestParams[which(row.names(object$forestParams) == "randomcombination"),i])
		if  ((randomCombinationString != "FALSE") )
		{	
			if (i == 1)
			{
				random.combination = as.numeric(unlist(strsplit(randomCombinationString, ",")))
				nbCombination = length(random.combination)/3
				X <- randomCombination(NAfactor2matrix(X, toGrep = "anythingParticular"), 
					combination = random.combination[1:(length(random.combination) - nbCombination)], 
					weights = random.combination[(length(random.combination) - nbCombination + 1):length(random.combination)])
			}
			else {	cat("Only first random combination will be used in incremental random uniform forest.\n") 	}
		}
	}	
	r.threads = threads	
	predObject <- randomUniformForestCore.predict(object, X, rf.aggregate = TRUE, pr.imbalance = classcutoff, pr.threads = threads, 	
	pr.classwt = classwt, pr.parallelPackage = parallelpackage[1])		
	object = filter.forest(object)	
	if (type[1] == "response") 
	{  
		predictedObject = predObject$majority.vote
		if (!is.null(classes)) 
		{  
			predictedObject = as.factor(predictedObject) 
			levels(predictedObject) = classes[as.numeric(levels(predictedObject))]
		}
	}	
	if (type[1] == "truemajority") 
	{
		{	
			nbObs = nrow(predObject$all.votes)
			trueMajorityVotes = rep(0,nbObs)
            alpha1 = (1 - conf)/2
			alpha2 = conf + alpha1			
			for (i in 1:nbObs)
			{
				expectedClasses = unlist(lapply(predObject$votes.data, function(Z) Z[i,1]))	
				votingMembers = unlist(lapply(predObject$votes.data, function(Z) Z[i,2]))
				if (object$regression)
				{
					votes = cbind(expectedClasses, votingMembers)
					colnames(votes) = c("expectedClasses", "votingMembers")
					trueMajorityVotes[i] = sum( votes[,"expectedClasses"]*votes[,"votingMembers"])/sum(votes[,"votingMembers"])
				}
				else
				{
					outliers = c(quantile(votingMembers, alpha1), quantile(votingMembers, alpha2))
					idxWoOutliers = which(votingMembers <= outliers[1] | votingMembers >= outliers[2])
					votes = cbind(expectedClasses[-idxWoOutliers], votingMembers[-idxWoOutliers])
					colnames(votes) = c("expectedClasses", "votingMembers")
					trueMajorityVotes[i] = which.max(by(votes, votes[,"expectedClasses"], sum))
				}
			}
			if (object$regression)	{ predictedObject = trueMajorityVotes }
			else {	predictedObject = as.factor(trueMajorityVotes); levels(predictedObject) = classes }
		}
	}	
	if (type[1] == "confInt")	
	{  
		if (!object$regression)
		{ stop( "confidence interval can only be computed for regression") }
		else
		{	
			alpha1 = (1 - conf)/2
			alpha2 = conf + alpha1
			Q_1 = apply(predObject$all.votes, 1, function(Z) quantile(Z, alpha1))
			Q_2 = apply(predObject$all.votes, 1, function(Z) quantile(Z, alpha2))
			SD = apply(predObject$all.votes, 1, function(Z) sd(Z))
			predictedObject = data.frame(cbind(predObject$majority.vote, Q_1, Q_2, SD))
			Q1Name = paste("LowerBound", "(q = ", round(alpha1,3), ")", sep ="")
			Q2Name = paste("UpperBound", "(q = ", round(alpha2,3), ")",sep ="")
			colnames(predictedObject) = c("Estimate", Q1Name, Q2Name, "standard deviation")
		}		
	}
	if (type[1] == "quantile")	
	{  
		if (!object$regression)
		{ stop( "quantile(s) can only be computed for regression") }
		else
		{	
			if (!is.numeric(whichQuantile)) 
			{ 
				stop( "Please use option 'whichQuantile', providing as argument a numeric value greater than 0 and lower than 1.\n") 
			}			
			if ( (whichQuantile <= 0) | (whichQuantile >= 1))  
			{ stop( "Please provide for 'whichQuantile' option, a numeric value than 0 and lower than 1.\n") }			
			predictedObject = apply(predObject$all.votes, 1, function(Z) quantile(Z, whichQuantile))
		}		
	}	
	if (type[1] == "all")	
	{  
		predictedObject = predObject
		if (object$regression)
		{ 
			stdDev = apply(predObject$all.votes, 1, sd)
			confIntObject = data.frame(cbind(predObject$majority.vote, apply(predObject$all.votes, 1, function(Z) quantile(Z,0.025)), 
			apply(predObject$all.votes, 1, function(Z) quantile(Z,0.975)), stdDev))			
			colnames(confIntObject) = c("Estimate", "LowerBound", "UpperBound", "Standard deviation")
			predictedObject$confidenceInterval = confIntObject
		}
	}	
	if (type[1] == "votes")  {	predictedObject = predObject$all.votes }	
	if ( (type[1] == "prob") & (!object$regression) )
	{	
		predictedObject = round(getVotesProbability2(predObject$all.votes, 1:length(classes)), 4) 
		colnames(predictedObject) = classes
	}	
	if ( (type[1] == "prob") & (object$regression) ) { stop("Probabilities can not be computed for regression") }	
	if ( (type[1] == "ranking") & (!object$regression) )
	{
		predictedObject = predictedObject2 = predObject$majority.vote
		if (!is.null(classes)) 
		{  
			if (!is.numeric(as.numeric(classes)))
			{	stop("Class are not numeric values or factor of numeric values. For Ranking, numeric values are needed as an equivalent of each class.") }
			else
			{
				minClass = min(as.numeric(classes))
				minPred = min(predictedObject)
				maxClass = max(as.numeric(classes))
				maxPred = max(predictedObject)				
				if (minClass < (minPred - 1))
				{	
					predictedObject[predictedObject != minPred] = predictedObject[predictedObject != minPred] - 1
					predictedObject = predictedObject - 1
				}
				else
				{
					if (minClass < minPred)
					{ 
				       predictedObject = predictedObject - 1
					}
				}		
				if (maxClass > maxPred)
				{	predictedObject[predictedObject == maxPred] = maxPred  }
			}
		}
		else { stop("Class not found. Please check if model is computed as classification, by looking if Y is set as factor.") }			
		numClasses = sort(unique(predObject$all.votes[,1]))
		probabilities = round(getVotesProbability2(predObject$all.votes, numClasses), 2) 
		colnames(probabilities) = classes
		countPredObject = table(predictedObject2)		
		majorityVote = as.numeric(names(which.max(countPredObject)))				
		n = length(predictedObject)
		followIdx = 1:n		
		if (!is.null(rankingIDs))
		{	
			rankingObject = cbind(followIdx, rankingIDs, predictedObject, probabilities)	
			if (length(classes) > 2)
			{  
				minoritiesProbability = rowSums(probabilities[,-majorityVote])
				rankingObject = cbind(rankingObject, minoritiesProbability)
			}			
			colnames(rankingObject)[1] = "idx"
			if (is.vector(rankingIDs) | is.factor(rankingIDs)) 
			{ 
				lengthRankingIDS = 1 
				colnames(rankingObject)[2] = "ID"
				colnames(rankingObject)[3] = "majority vote"
				cases = sort(unique(rankingIDs))
			}
			else 
			{ 
				lengthRankingIDS = ncol(rankingIDs)
				if (is.null(colnames(rankingIDs)))
				{	colnames(rankingObject)[2:(2+lengthRankingIDS)] = "ID"	}			
				
				colnames(rankingObject)[lengthRankingIDS+2] = "majority vote"
				cases = sort(unique(rankingIDs[,1]))
			}			
			if (length(classes) > 2)
			{  minorityIdx = which(colnames(rankingObject) == "minoritiesProbability") }
			else
			{ 	minorityIdx = which(colnames(rankingObject) == classes[-majorityVote])	}			
			lengthCases = length(cases)
			subCases = vector('list', lengthCases)
			for (i in 1:lengthCases)
			{	subCases[[i]] = which(rankingObject[,2] == cases[i]) }							
			rankingOutputObject <- matrix(NA, n, ncol(rankingObject) + 1)
			for (i in 1:lengthCases)
			{
				if (length(subCases[[i]]) > 1)
				{
					rankingOutputObject[subCases[[i]],] <- as.matrix(cbind(sortMatrix(rankingObject[subCases[[i]],], minorityIdx, decrease = TRUE), 
					1:length(subCases[[i]])))
				}
				else
				{	rankingOutputObject[subCases[[i]],] <- c(rankingObject[subCases[[i]],], 1)	}
			}			
			rankingOutputObject = as.data.frame(rankingOutputObject)			
			for (j in 1:ncol(rankingOutputObject))
			{
				if (is.factor(rankingOutputObject[,j]) & is.numeric(as.numeric(as.vector(rankingOutputObject[,j]))))
				{	rankingOutputObject[,j] = as.numeric(as.vector(rankingOutputObject[,j])) }
			}			
			colnames(rankingOutputObject) = c(colnames(rankingObject), "rank")
			predictedObject = sortMatrix(rankingOutputObject,1)	
		}
		else
		{	
			rankingObject = cbind(followIdx, predictedObject, probabilities)			
			if (length(classes) > 2)
			{  
				minoritiesProbability = rowSums(probabilities[,-majorityVote])
				rankingObject = cbind(rankingObject, minoritiesProbability)
			}			
			colnames(rankingObject)[1] = "idx"
			colnames(rankingObject)[2] = "majority vote"			
			if (length(classes) > 2)
			{  minorityIdx = which(colnames(rankingObject) == "minoritiesProbability") }
			else
			{	minorityIdx = which(colnames(rankingObject) == classes[-majorityVote]) }			
			rankingOutputObject = sortMatrix(rankingObject, minorityIdx, decrease = TRUE)	
			predictedObject = cbind(rankingOutputObject, 1:n)
			colnames(predictedObject)[ncol(predictedObject)] = "rank"
			predictedObject = sortMatrix(predictedObject,1)	
		}
	}		
	if ( (type[1] == "ranking") & (object$regression) )
	{ stop("Ranking currently available for classification tasks") }	
	predictedObject
}

randomUniformForestCore.predict <- function(object, X, 
rf.aggregate = TRUE, 
OOB.idx = FALSE, 
pr.imbalance = c(0,0),
pr.regression = TRUE,
pr.classwt = FALSE,  
classes = -1,
pr.export = c("onlineClassify", "predictDecisionTree", "majorityClass", "randomWhichMax", "rmNA", "mergeLists", "classifyMatrixCPP", "L2DistCPP", "checkUniqueObsCPP", "crossEntropyCPP", "giniCPP", "L2InformationGainCPP",	"entropyInformationGainCPP", "runifMatrixCPP", "NATreatment"), 
pr.threads = "auto", 
pr.parallelPackage = "doParallel")
{
	## #require(rUniformForestCppClass)	
	if (!OOB.idx)
	{
		if (is.matrix(X))
		{
			matchNA = (length(which(is.na(X))) > 0)		
			if (matchNA) 
			{ 
				cat("NA found in data. Fast imputation (means) is used for missing values. Please use one of many models available if accuracy is needed\n Or use na.impute() function with the option 'na.action = accurateImpute'.")
				X <- na.impute(X)
			}
		}
		else
		{ 
			X <- NAfactor2matrix(X, toGrep = "anythingParticular")	
			matchNA = (length(which(is.na(X))) > 0)		
			if (matchNA) 
			{ 
				cat("NA found in data. Fast imputation (means) is used for missing values. Please use one of many models available if accuracy is needed\n Or use na.impute() function with the option 'na.action = accurateImpute'.\n")
				X <- na.impute(X)
			}
		}
		if (!is.null(object$logX))
		{
			if (object$logX)
			{ 
			   if (is.null(object$categoricalvariables))  {  X <- generic.log(X) }
			   else  {  X[,-object$categoricalvariables] <- generic.log(X[,-object$categoricalvariables]) }
			}
		}		
		object = filter.forest(object)
	}				
	if (!is.null(object$OOB) & (!OOB.idx))
	{  OOB.predicts = object$OOB.predicts; pr.regression = object$regression; object = object$object }
	else
	{
		if (!OOB.idx) 	{ 	pr.regression = object$regression; object = object$object   }			
		if (!is.null(object$variableImportance) & (OOB.idx)) {	object = object$object	}		
	}			
	n = nrow(X)			
	if (!pr.regression)
	{
		if (classes[1] <  1)
		{
			classes = sort(as.numeric(rownames(table(object[[sample(1:length(object),1)]][,"prediction"]))))
			if (classes[1] == 0) {  classes = classes[-1]  }
		}			
		l.class = length(classes)
		class.occur = rep(0,l.class)
	}	
	{
		ntree = length(object)
		pred.X = vector("list", ntree)
		all.votes = nodes.length = nodes.depth = matrix(data = Inf, nrow = n, ncol = ntree)
		majority.vote = vector(length = n)		
		fullDim = ntree*n		
		if ((fullDim < 1e7) | (pr.threads == 1))
		{	pred.X <- lapply(object, function(Z) predictDecisionTree(Z, X))	 }
		else
		{
			{
				## #require(parallel)	
				max_threads = detectCores()				
				if (pr.threads == "auto")
				{	
					if (max_threads == 2) { pr.threads = max_threads }
					else {	pr.threads  = max(1, max_threads - 1)  }
				}
				else
				{
					if (max_threads < pr.threads) 
					{	cat("Note: number of threads is higher than logical threads in this computer.\n") }
				}
				## #require(doParallel)	
				Cl = makePSOCKcluster(pr.threads, type = "SOCK")
				registerDoParallel(Cl)
				chunkSize  <-  ceiling(n/getDoParWorkers())
				smpopts  <- list(chunkSize = chunkSize)
			}						
			pushObject <- function(object, X) lapply(object, function(Z) predictDecisionTree(Z, X))			
			pred.X <- foreach(X = iterators::iter(X, by ='row', chunksize = chunkSize), .combine = mergeLists, .export = pr.export) %dopar%
			pushObject(object, X)									
			stopCluster(Cl)
		}								
		for (i in 1:ntree)
		{	
			all.votes[,i] = pred.X[[i]][,1]
			nodes.length[,i] = pred.X[[i]][,2]
			nodes.depth[,i] = pred.X[[i]][,3]	
		}		
		if (pr.classwt)
		{
			allWeightedVotes = matrix(data = Inf, nrow = n, ncol = ntree)
			for (i in 1:ntree)	{  allWeightedVotes[,i] = object[[i]][nodes.depth[,i], "avgLeafWeight"]	}
		}
		else
		{  allWeightedVotes = NULL }
			
		if (pr.regression) {  majority.vote = rowMeans(all.votes)	}
		else {	majority.vote = majorityClass(all.votes, classes, m.imbalance = pr.imbalance, m.classwt = allWeightedVotes)$majority.vote  }
	}				
	if (rf.aggregate & (!OOB.idx))
	{	
		return(list(all.votes = all.votes, majority.vote = majority.vote, nodes.depth = nodes.depth, nodes.length = nodes.length, 
		votes.data = pred.X))	
	}	
	else
	{
		if (OOB.idx) 
		{  
			if (pr.classwt)	{ return(list(all.votes = all.votes, allWeightedVotes = allWeightedVotes))	}
			else { return(all.votes)	}
		}
		else  {  return(majority.vote) 	}		
	}
}

majorityClass <- function(all.votes, classes, m.regression = FALSE, m.imbalance = c(0,0), m.classwt = NULL, 
m.threads = "auto")
{
	if (m.regression)
	{	
		all.votes[all.votes == Inf] = NA  
		majority.vote = rowMeans(all.votes, na.rm = TRUE) 
		class.counts = NULL	
		return(list(majority.vote = majority.vote, class.counts = class.counts))
	}
	else
	{
		n = nrow(all.votes)
		majority.vote = vector(length = n)
		l.class = length(classes)
		class.occur = vector(length = l.class)
		class.countsMajorityVote = trueClass.counts = matrix(NA, n, l.class)
		for (i in 1:n)
		{ 
			if (!is.null(m.classwt))
			{	
				for (j in 1:l.class) 
				{
					idx = rmNA(which(all.votes[i,] == classes[j]))
					l.idx = length(idx)									
					if (l.idx > 0) {   class.occur[j] = sum(m.classwt[i,idx]) }
					else	{ class.occur[j] = 0 }
				}				
			}
			else
			{
				for (j in 1:l.class) { 	class.occur[j] <- sum(all.votes[i,] == classes[j])  }
			}				
			if (m.imbalance[1] != 0)
			{	class.occur[m.imbalance[1]] = floor(class.occur[m.imbalance[1]]*m.imbalance[2]) } 
			majority.vote[i] = randomWhichMax(class.occur)
			class.countsMajorityVote[i,] = c(class.occur[majority.vote[i]], class.occur[-majority.vote[i]])
			trueClass.counts[i,]  = class.occur
		}
		return(list(majority.vote = majority.vote, class.counts = class.countsMajorityVote, trueClass.counts = trueClass.counts))
	}	
}

getVotesProbability  <- function(X, classes) (majorityClass(X, classes)$class.counts/ncol(X))

getVotesProbability2  <- function(X, classes) (majorityClass(X, classes)$trueClass.counts/ncol(X))

strength_and_correlation <- function(OOB.votes, OOB.object, 
	rf.classes, 
	s.trainLabels = NULL, 
	s.regression = FALSE, 
	output = NULL, 
	s.threads = "auto", 
	s.parallelPackage = "doParallel" )
{
	j = NULL
	dimOOBvotes = dim(OOB.votes)	
	n.OOB = dimOOBvotes[1]
	p.OOB = dimOOBvotes[2]	
	if (s.regression)
	{
		OOB.votes[OOB.votes == Inf] = NA		
		if (is.null(s.trainLabels))  { Y = rowMeans(OOB.votes, na.rm =TRUE)	}
		else {  Y = s.trainLabels  }			
		expectedSquaredErrorOverTrees = colMeans(  (OOB.votes - Y)^2, na.rm = TRUE)
		PE.tree = sum(expectedSquaredErrorOverTrees)/length(expectedSquaredErrorOverTrees)
		sd.T =  sqrt(expectedSquaredErrorOverTrees)
		mean.corr = mean(rowMeans(cor((Y - OOB.votes), use = "pairwise.complete.obs")))		
		if (is.na(mean.corr))
		{	
			cat("Not enough data to compute average correlation for trees. Error is then prediction error of a tree.\n")
			return(list(PE.forest = (mean(sd.T))^2,  PE.max = PE.tree, PE.tree = PE.tree, mean.corr = mean.corr))
		}
		else
		{	return(list(PE.forest = mean.corr*(mean(sd.T))^2,  PE.max = mean.corr*PE.tree, PE.tree = PE.tree, mean.corr = mean.corr))	}
	}
	else
	{
		## #require(parallel)		
		max_threads = detectCores()		
		if (s.threads == "auto")
		{	s.threads  = max(1, max_threads - 1)  }
		else
		{
			if (max_threads < s.threads) 
			{	cat("Note: number of threads is higher than logical threads in this computer.\n") }
		}
		{	
			## #require(doParallel)
			Cl <- makePSOCKcluster(s.threads, type ="SOCK")
			registerDoParallel(Cl)
		}
		chunkSize  <-  ceiling(p.OOB/getDoParWorkers())
		smpopts  <-  list(chunkSize = chunkSize)	
		p.new.OOB = apply(OOB.votes, 1, function(OOB.votes) sum(OOB.votes != Inf))			
		if (length(rf.classes) == 2)
		{
			Q.x.1 = OOB.object$class.counts[,1]/rowSums(OOB.object$class.counts)
			rawStrength = 2*Q.x.1 - 1			
			Tk.1 = Tk.2 = matrix (data = NA, ncol = p.OOB, nrow = n.OOB)
			OOB.votes.1 = cbind(OOB.votes, OOB.object$majority.vote)
			Tk.1 <- foreach(j = 1:p.OOB, .options.smp = smpopts, .combine = cbind, .multicombine = TRUE) %dopar%
			apply(OOB.votes.1[,-j], 1, function(Z) sum(Z[-p.OOB] == Z[p.OOB]))
			Tk.1 = Tk.1/p.new.OOB 
			Tk.2 = 1 - Tk.1
			stopCluster(Cl)
		}
		else
		{
			nn = rowSums(OOB.object$class.counts)
			Q.x.j = apply(OOB.object$class.counts[,-1], 1, max)
			Q.x.y = OOB.object$class.counts[,1]			
			rawStrength = (Q.x.y - Q.x.j)/nn
			maj.class.j = vector(length = n.OOB)
			for (i in 1:n.OOB)
			{  
				second.class.max = randomWhichMax(OOB.object$class.counts[i,-1])
				maj.class.j[i] = rf.classes[-OOB.object$majority.vote[i]][second.class.max]
			}
			OOB.votes.1 = cbind(OOB.votes, OOB.object$majority.vote)
			OOB.votes.2 = cbind(OOB.votes,  maj.class.j)			
			ZZ <- function(j) cbind(apply(OOB.votes.1[,-j], 1, function(Z) sum(Z[-p.OOB] == Z[p.OOB])), apply(OOB.votes.2[,-j], 1, 
				function(Z) sum(Z[-p.OOB] == Z[p.OOB])))					
			Tk <- foreach(j = 1:p.OOB, .options.smp = smpopts, .combine = cbind, .multicombine = TRUE) %dopar% ZZ(j)		
			mixIdx = getOddEven(1:ncol(Tk))			
			Tk.1 = Tk[,mixIdx$odd]
			Tk.2 = Tk[,mixIdx$even]
			Tk.1 = Tk.1/p.new.OOB
			Tk.2 = Tk.2/p.new.OOB
			stopCluster(Cl) 
		}				
		p1 = colMeans(Tk.1)
		p2 = colMeans(Tk.2)		
		strength = mean(rawStrength)
		varStrength =  var(rawStrength)
		sd.T = ((p1 + p2 + (p1 - p2)^2))^0.5		
		mean.corr = varStrength / (mean(sd.T))^2			
		PE.est = mean.corr*(1 -	strength^2)/strength^2		
		return(list(PE = PE.est, avg.corr = mean.corr, strength = strength, std.strength = sqrt(varStrength)))
	}
}

monitorOOBError  <- function(OOB.votes, Y, regression = FALSE, threads = "auto", f = L2Dist)
{
	j = NULL
	n = length(Y)
	n.RS = sqrt(n)
	p = ncol(OOB.votes)	
	## #require(parallel)	
	max_threads = detectCores()		
	if (threads == "auto")
	{	
		if (max_threads == 2) { threads = max_threads }
		else {	threads  = max(1, max_threads - 1)  }
	}
	else
	{
		if ((max_threads) < threads) 
		{	 cat("Note: number of threads is higher than logical threads in this computer\n.") }
	}	
	## #require(doParallel)			
	Cl = makePSOCKcluster(threads, type = "SOCK")
	registerDoParallel(Cl)	
	chunkSize  <-  ceiling(p/getDoParWorkers())
	smpopts  <- list(chunkSize  =  chunkSize)
	if (!regression)
	{
		Y = as.numeric(Y)
		classes = sort(unique(Y))
		S1 = sample(classes,1)
		S2 = classes[-S1][1]
		
		OOBmonitor <- foreach(j = 1:(p-1), 
		.export = c("generalization.error", "confusion.matrix", "majorityClass", "rmNA", "randomWhichMax"), 
		.options.smp = smpopts, .combine = rbind) %dopar%
		{
			Estimate <- generalization.error(confusion.matrix(majorityClass(OOB.votes[,1:(j+1)], classes)$majority.vote, Y))
			C1 <- generalization.error(confusion.matrix(majorityClass(OOB.votes[,1:(j+1)], classes, m.imbalance =c(1, 1.5))$majority.vote, Y))
			C2 <- generalization.error(confusion.matrix(majorityClass(OOB.votes[,1:(j+1)], classes, m.imbalance= c(2, 1.5))$majority.vote, Y))
			
			t(c(Estimate, C1, C2))
		}
		
		E0 <- generalization.error(confusion.matrix(majorityClass(matrix(OOB.votes[,1]), classes)$majority.vote, Y))
		C1 <- generalization.error(confusion.matrix(majorityClass(matrix(OOB.votes[,1]), classes, m.imbalance =c(1, 1.5))$majority.vote, Y))
		C2 <- generalization.error(confusion.matrix(majorityClass(matrix(OOB.votes[,1]), classes, m.imbalance= c(2, 1.5))$majority.vote, Y))
		OOBmonitor = rbind(t(c(E0, C1, C2)),OOBmonitor)
	}
	else
	{
		OOBmonitor <- foreach(j = 1:(p-1), .export = c("generalization.error", "confusion.matrix", "majorityClass", "rmNA", "L2Dist", "L1Dist"), .options.smp = smpopts, .combine = c) %dopar%
		{
			Z = majorityClass(OOB.votes[,1:(j+1)], 0, m.regression = TRUE)$majority.vote
			NAIdx = which(is.na(Z))
			if (length(NAIdx) > 0) { f(Z[-NAIdx], Y[-NAIdx])/length(Z[-NAIdx])}
			else { f(Z, Y)/length(Z)  }
		}		
		Z = majorityClass(matrix(OOB.votes[,1]), 0, m.regression = TRUE)$majority.vote
		NAIdx = which(is.na(Z))
		E0 = if (length(NAIdx) > 0) { f(Z[-NAIdx], Y[-NAIdx])/length(Z[-NAIdx])} else { f(Z, Y)/length(Z)  }
		OOBmonitor = c(E0,OOBmonitor)
	}	
	stopCluster(Cl)
	return(OOBmonitor)
}

weightedVote <- function(all.votes, idx = 2, granularity = 2)
{
	all.votes <- round(all.votes, granularity)	
	apply(all.votes, 1, function(Z)
		{
			A = sort(table(rmInf(Z)), decreasing = TRUE)
			B = as.numeric(names(A))[1:idx]
			sum( B * (B/sum(abs(B))) )
		}
	)
}

weightedVoteModel <- function(votes, majorityVote, Y = NULL, nbModels = 1, idx = 1 , granularity = 1, train = TRUE, models.coeff = NULL)
{
	if (train)
	{
		if (is.null(Y))
		{ stop("output is neded to build model") }		
		if (nbModels == 1)
		{
			default.model = weightedVote(votes, idx = idx, granularity = granularity)
			lin.model = lm(Y ~ cbind(majorityVote, default.model))
		}
		else
		{
			models = matrix(data = NA, ncol = nbModels + 2, nrow = length(Y))			
			models[,1] = majorityVote
			models[,2] = weightedVote(votes, idx = idx, granularity = granularity)  
			for (j in 3:(nbModels+2))
			{	models[,j] = weightedVote(votes, idx = max(j,idx), granularity = j ) }
			lin.model = lm(Y ~ models)
		}		
		return(lin.model)
	}
	else
	{
		p = length(models.coeff)
		models = matrix(data = NA, ncol = p, nrow = length(majorityVote))			
		models[,1] = rep(1,length(majorityVote))
		models[,2] = majorityVote
		models[,3] = weightedVote(votes, idx = idx, granularity = granularity)		
		if (p > 3)
		{
			for (j in 4:p)
			{	models[,j] = weightedVote(votes,  idx =j  , granularity = j)	}
		}		
		newMajorityVote = apply(models,1, function(Z) sum(models.coeff*Z))		
		return(newMajorityVote )
	}
}

postProcessingVotes <- function(object, nbModels = 1, idx = 1, granularity = 1, predObject = NULL, swapPredictions = FALSE, X = NULL, imbalanced = FALSE)
{
	object <- filter.object(object)	
	if (rownames(object$forestParams)[1] == "reduceDimension") 
	{ stop("Post Processing does not work with objects coming from rUniformForest.big() function") }
	if (!object$forest$regression)
	{
		if (length(as.numeric(object$classes)) > 2)
		{	stop("Optimization currently works only for binary classification") }		
		if (is.null(X))
		{	stop("Please provide test data")	}
		else
		{
			if (is.null(predObject)) {	predObject <- predict(object, X, type = "all")	}
			else
			{
				if (is.null(predObject$all.votes)) 
				{ stop("Please provide full prediction object (option type = 'all' when calling predict() function).") }				
				if (is.null(predObject$majority.vote)) 
				{ stop("Please provide full prediction object (option type = 'all' when calling predict() function).") }
			}
			majorityVotePosition = which.max(table(predObject$majority.vote))
			numClasses = sort(unique(predObject$all.votes[,1]))
			probPred = round(getVotesProbability2(predObject$all.votes, numClasses), 2) 
			colnames(probPred) = object$classes
			if (imbalanced) 
			{ 	cutoff =  1 - ( mean(probPred[,2])/mean(probPred[,1]) )	}
			else 
			{ 	cutoff = 0.5/2*( mean(probPred[,2])/mean(probPred[,1]) + mean(probPred[,1])/mean(probPred[,2]) ) } 
			predVotes = predict(object, X, classcutoff = c(object$classes[majorityVotePosition], cutoff))			
			return(predVotes)
		}
	}	
	if (swapPredictions) 
	{ 
		object$forest$OOB.votes = predObject$forest$OOB.votes 
		object$forest$OOB.predicts = predObject$forest$OOB.predicts
	}	
	if (is.null(object$forest$OOB.votes)) 
	{ stop("No OOB data for post processing. Please enable OOB option and subsamplerate or bootstrap ('replace' option) when computing model.") }			
	if (is.null(object$predictionObject))
	{
		if (is.null(predObject)) { stop("Post processing can not be computed. Please provide prediction object") }
		else 
		{  
			if (is.null(predObject$majority.vote))
			{ stop("Post processing can not be computed. Please provide full prediction object (type = 'all') when calling predict()") }
			else
			{			
				meanEstimate = predObject$majority.vote
				allVotes = predObject$all.votes
			}
		}
	}
	else
	{   
		meanEstimate = object$predictionObject$majority.vote
		allVotes = object$predictionObject$all.votes
	}			
	Y = object$y
	OOBVotes = object$forest$OOB.votes	
	NAIdx = which.is.na(object$forest$OOB.predicts)	
	if (length(NAIdx) > 0)
	{
		Y = Y[-NAIdx]
		OOBVotes = OOBVotes[-NAIdx,]
		object$forest$OOB.predicts = object$forest$OOB.predicts[-NAIdx]
	}
	OOBMeanEstimate = object$forest$OOB.predicts
	L2DistOOBMeanEstimate <- L2Dist(OOBMeanEstimate, Y)/length(Y)	
	OOBMedianEstimate <- apply(OOBVotes, 1, function(Z) median(rmInf(Z)))
	L2DistOOBMedianEstimate <- L2Dist(OOBMedianEstimate, Y)/length(Y)	
	weightedModel <- weightedVoteModel(OOBVotes, OOBMeanEstimate, Y = Y, nbModels = nbModels, idx = idx, granularity = granularity)
	OOBWeightedVoteEstimate <- weightedVoteModel(OOBVotes, OOBMeanEstimate, train = FALSE, models.coeff = weightedModel$coefficients)	
	NAOOBWeightedVoteEstimate <- which.is.na(OOBWeightedVoteEstimate)	
	if (length(NAOOBWeightedVoteEstimate) > 0) 
	{ OOBWeightedVoteEstimate[NAOOBWeightedVoteEstimate] = OOBMeanEstimate[NAOOBWeightedVoteEstimate] }	
	L2DistOOBWeightedVoteEstimate <- L2Dist(OOBWeightedVoteEstimate, Y)/length(Y)	
	flagMedian =  ( (L2DistOOBMedianEstimate <= L2DistOOBMeanEstimate) | (L1Dist(OOBMedianEstimate, Y) <= L1Dist(OOBMeanEstimate, Y)) )
	flagWeighted = ( (L2DistOOBWeightedVoteEstimate <= L2DistOOBMeanEstimate) | (L1Dist(OOBWeightedVoteEstimate, Y) <= L1Dist(OOBMeanEstimate, Y)) )	
	if (flagMedian)
	{
		if (flagWeighted)
		{  weightedVoteEstimate <- weightedVoteModel(allVotes, meanEstimate, train = FALSE, models.coeff = weightedModel$coefficients) }
		else
		{  weightedVoteEstimate = rep(0, length(meanEstimate)) }
	}
	else
	{  	weightedVoteEstimate <- weightedVoteModel(allVotes, meanEstimate, train = FALSE, models.coeff = weightedModel$coefficients)  }	
	medianEstimate <- apply(allVotes, 1, function(Z) median(rmInf(Z)))	
	NAWeightedVoteEstimate <- which.is.na(weightedVoteEstimate)
	if (length(NAWeightedVoteEstimate) > 0) { weightedVoteEstimate[NAWeightedVoteEstimate] = meanEstimate[NAWeightedVoteEstimate] }	
	negativeBias = - (mean(OOBMeanEstimate) - Y)
	modelResiduals = Y - OOBMeanEstimate
	biasModel  = lm(modelResiduals ~ negativeBias) 
	biasSign = sign(meanEstimate - medianEstimate)   
	biasSign[which(biasSign == 0)] = 1
	residuals_hat = vector(length = ncol(OOBVotes))	
	residuals_hat <- apply(OOBVotes, 2, function(Z) sum(rmInf(Z))/length(rmInf(Z))) - mean(OOBMeanEstimate)
	MeanNegativeBias = -mean(residuals_hat^2)		
	biasCorrection = biasModel$coefficients[2]*biasSign*MeanNegativeBias + biasModel$coefficients[1]	
	if (flagMedian)
	{ 
		if (flagWeighted) { return( 0.5*(medianEstimate + weightedVoteEstimate + biasCorrection)) }
		else { return(medianEstimate) }
	}
	else
	{ return(weightedVoteEstimate + biasCorrection) }
}

localVariableImportance <- function(object, nbVariables = 2, Xtest = NULL, predObject = NULL, l.threads = "auto", l.parallelPackage = "doParallel")
{
	object = filter.object(object)
	rufObject = object$forest$object	
	if (is.null(predObject)) 
	{ 
		if (is.null(Xtest))
		{ stop("Local variable importance can not be computed. please provide test data.") }
		else
		{
			predObject = predict(object, Xtest, type = "all")
			majority.vote = as.numeric(predObject$majority.vote) 
			pred.rufObject = predObject$votes.data
		}
	}
	else 
	{  
		if (is.null(predObject$majority.vote))
		{ 
			stop("Local variable importance can not be computed. Please provide full prediction object (type = 'all') when calling predict()") 
		}
		else
		{			
			majority.vote = as.numeric(predObject$majority.vote) 
			pred.rufObject = predObject$votes.data
		}		
	}
	ntree = length(rufObject)
	if (is.null(pred.rufObject))
	{ stop ("no data to evaluate importance") }
	else
	{
		n = dim(pred.rufObject[[1]])[1]
		{
			## #require(parallel)
			threads = l.threads
			max_threads = detectCores()			
			if (threads == "auto")
			{	
				if (max_threads == 2) { threads = max_threads }
				else {	threads  = max(1, max_threads - 1)  }
			}
			else
			{
				if ((max_threads) < threads) 
				{	cat("Note: number of threads is higher than logical threads in this computer.\n") }
			}
			{
				## #require(doParallel)				
				Cl = makePSOCKcluster(threads, type = "SOCK")
				registerDoParallel(Cl)
			}
			chunkSize  <- ceiling(ntree/getDoParWorkers())
			smpopts  <- list(chunkSize = chunkSize)
		}
		b = NULL
		varMatrix <- foreach(b = 1:ntree, .options.smp = smpopts, .combine = rbind, .multicombine = TRUE) %dopar%	
		{ 	
			predVar = vector(length = n)
			for (i in 1:n)
			{
				predIdx = pred.rufObject[[b]][i,3]
				predVar[i] = rufObject[[b]][which(rufObject[[b]][,"left daughter"] == predIdx 
				| rufObject[[b]][,"right daughter"] == predIdx), "split var"][1]
			}						
			predVar
		}
		stopCluster(Cl)
		varImp = varFreq = matrix(data = NA, ncol = nbVariables, nrow = n)
		varImp <- t(apply(varMatrix, 2, function(Z) 
			as.numeric(names(sort(table(Z), decreasing = TRUE)[1:nbVariables])) ))
		varFreq <- t(apply(varMatrix, 2, function(Z)  sort(table(Z), decreasing = TRUE)[1:nbVariables]))/ntree	
		NAIdx = which(is.na(varImp), arr.ind = TRUE)		
		if (length(NAIdx[,1]) > 0)
		{	
			varImp[NAIdx] = apply(NAIdx, 1, function(Z) { 
												newZ = rmNA(varImp[Z])
												if (length(newZ) == 1) { rep(newZ, length(Z)) } else  { sample(newZ,1)  }
											}
							)
			varFreq[NAIdx] = apply(NAIdx, 1, function(Z) { 
												newZ = rmNA(varFreq[Z])
												if (length(newZ) == 1) { rep(newZ, length(Z)) } else  { sample(newZ,1)  }
											}
							)				
		}		
		objectNames = objectFrequencyNames = vector()
		for (j in 1:nbVariables)
		{	
			objectNames[j] = paste("localVariable", j, sep = "")
			objectFrequencyNames[j] = paste("localVariableFrequency", j, sep = "") 
		}	
		variableImportance.object = cbind(majority.vote, varImp, varFreq)				
		if (!object$forest$regression)
		{
			colnames(variableImportance.object) = c("class", objectNames, objectFrequencyNames)		
			classObject = as.numeric(names(table(variableImportance.object[,1])))		
			classMatrix = list()		
			for (i in 1:length(classObject))
			{		
				classMatrix[[i]] = round(sort(table(variableImportance.object[which(variableImportance.object[,1] == i),2])/
					sum(table(variableImportance.object[which(variableImportance.object[,1] == i),2])), decreasing = TRUE),2)
			}		
			orderedVar = unique(as.numeric(names(sort(unlist(classMatrix), decreasing = TRUE))))
			Variables = as.numeric(names(sort(unlist(classMatrix), decreasing = TRUE)))
			orderedImportance =  as.numeric(sort(unlist(classMatrix), decreasing = TRUE))			
			classVariableImportance = matrix(data = 0, ncol = length(classObject), nrow = length(orderedVar))
			for (i in 1:length(orderedVar))
			{	
				for (j in 1:length(classObject))
				{	
					idx = which( as.numeric(names(classMatrix[[j]])) == orderedVar[i])
					if (length(idx) > 0)
					{ classVariableImportance[i,j] = classMatrix[[j]][idx]	}
				}
			}		
			classVariableImportance = data.frame(classVariableImportance)
			rownames(classVariableImportance) = orderedVar
			colnames(classVariableImportance) = paste("Class", classObject, sep=" ")				
			return(list( obsVariableImportance = variableImportance.object, classVariableImportance = classVariableImportance ))
		}
		else
		{
			colnames(variableImportance.object) = c("estimate", objectNames, objectFrequencyNames)
			return(variableImportance.object)			
		}
	}
}
	
localTreeImportance <- function(rfObject, OOBPredicts = NULL, OOBVotes = NULL)
{
	rfObject = filter.object(rfObject)  
	if (!is.null(rfObject$OOB.votes))	{ OOBVotes = rfObject$OOB.votes }
	if (is.null( OOBVotes))	{ stop ( "OOB outputs missing. No optimization or weighted vote can be done") }
	if (!is.null(rfObject$OOB.predicts))	{	OOBPredicts = rfObject$OOB.predicts	}	
	n = length(OOBPredicts)  
	wPlus = wMinus = NULL
	for ( i in seq_along(OOBPredicts))
	{
		wPlus = c(wPlus, which(OOBVotes[i,] == OOBPredicts[i]))
		wMinus = c(wMinus, which(OOBVotes[i,] != OOBPredicts[i]))
	}	
	object = list(pweights = (table(wPlus)/n), nweights = (table(wMinus)/n))	
	return(object)
}

importance.randomUniformForest <- function(object, maxVar = 30, maxInteractions = 3, Xtest = NULL, 
predObject = NULL, ...)
{
	object <- filter.object(object) 
	maxInteractions = max(2, maxInteractions)	
	if (!is.null(Xtest))
	{
		if (!is.null(object$formula) & (length(object$variablesNames) != ncol(Xtest)))
		{
			mf <- model.frame(formula = object$formula, data = as.data.frame(Xtest))
			Xtest <- model.matrix(attr(mf, "terms"), data = mf)[,-1]
			if (object$logX) { 	Xtest <- generic.log(Xtest) }		
		}
		else
		{		
			Xfactors <- which.is.factor(Xtest)
			catVarIdx <- which(rownames(object$paramsObject) == "categorical variables")
			Xtest <- NAfactor2matrix(Xtest, toGrep = "anythingParticular")
			matchNA = (length(which(is.na(Xtest))) > 0)				
			if (matchNA) 
			{ 
				cat("NA found in data. Fast imputation (means) is used for missing values\n")
				Xtest <- na.impute(Xtest)
			}			
			if (object$logX)
			{ 	
				trueFactorsIdx <- which(Xfactors == 1)
				if (length(trueFactorsIdx) ==  0)  {  Xtest <- generic.log(Xtest) }
				else  
				{  
					if (!is.null(object$paramsObject[1,catVarIdx]))
					{	
						cat("Warning : all categorical variables have been ignored for the logarithm transformation.\nIf only some of them have been defined as categorical, variable importance will be altered.\n To avoid it, please define logarithm transformation outside of the model or ignore categorical variables, or set them to 'all'.\n")
						Xtest[,-trueFactorsIdx] <- generic.log(Xtest[,-trueFactorsIdx])						
					}
					else
					{ 	Xtest <- generic.log(Xtest)	}
				}
			}
		}		
	}	
	if (!is.null(object$forest$variableImportance))
	{	
		par(las=1) 
		maxChar = floor(1 + max(nchar(object$variablesNames))/2)		
		par(mar=c(5, maxChar + 1,4,2)) 		
		varImportance1 = varImportance = object$forest$variableImportance
		if (!object$forest$regression)
		{
			varImportance[,"class"] = object$classes[as.numeric(varImportance[,"class"])]
			varImportance1[,"class"] = varImportance[,"class"]
		}		
		nVar = nrow(varImportance)
		if (nVar > maxVar) { varImportance = varImportance[1:maxVar,] }	
		barplot(varImportance[nVar:1, "percent.importance"], horiz = TRUE, col = sort(heat.colors(nVar), 
		decreasing = TRUE), 
		names.arg = varImportance[nVar:1,"variables"], xlab = "Relative Influence (%)", 
		main = if (object$forest$regression)  { "Variable Importance based on 'Lp' distance" }  
		else  { "Variable Importance based on information gain" }, border = NA)
		abline(v = 100/nVar, col = 'grey')				
		localVarImportance <- localVariableImportance(object, nbVariables = maxInteractions, Xtest = Xtest, 
		predObject = predObject)		
		if (!object$forest$regression)
		{
			cat("\n1 - Global Variable Importance (", min(maxVar, nVar), " most important based on information gain) :\n", sep = "")
			cat("Note: most predictive features are ordered by 'score' and plotted. Most discriminant ones\nshould also be taken into account by looking 'class' and 'class.frequency'.\n\n")				
			rownames(localVarImportance$classVariableImportance) = object$variablesNames[as.numeric(rownames(localVarImportance$classVariableImportance))]
			colnames(localVarImportance$classVariableImportance) = paste("Class ", 
				object$classes[as.numeric(rm.string(colnames(localVarImportance$classVariableImportance), "Class"))], sep = "")				
			obsVarImportance = data.frame(localVarImportance$obsVariableImportance)
			obsVarImportance[,1] = object$classes[obsVarImportance[,1]]
		}
		else
		{	
			cat("\n1 - Global Variable Importance (", min(maxVar, nVar), 
			" most important based on 'Lp' distance) :\n", sep = "")
			obsVarImportance = data.frame(localVarImportance) 	
		}
		print(varImportance)		
		for (j in 2:(maxInteractions+1)) { obsVarImportance[,j] = object$variablesNames[obsVarImportance[,j]] }		
		obsVarImportance2 = obsVarImportance
		fOrder = sort(table(obsVarImportance2[,2]), decreasing = TRUE)
        sOrder = sort(table(obsVarImportance2[,3]), decreasing = TRUE)				
		W1 = mean(obsVarImportance2[,grep("localVariableFrequency1", colnames(obsVarImportance2))[1]])
		W2 = mean(obsVarImportance2[,grep("localVariableFrequency2", colnames(obsVarImportance2))[1]])		
		minDim2 = min(length(fOrder), length(sOrder))
		partialDependence = matrix(NA, minDim2, minDim2)
		for (i in 1:minDim2) {  partialDependence[,i] = fOrder[i]*W1/W2 + sOrder[1:minDim2]  }		
		colnames(partialDependence) = names(fOrder)[1:minDim2]
		rownames(partialDependence) = names(sOrder)[1:minDim2]
		partialDependence = partialDependence/(2*nrow(obsVarImportance2))
		avg1rstOrder = colMeans(partialDependence)
		avg2ndOrder = c(rowMeans(partialDependence),0)
		partialDependence = rbind(partialDependence, avg1rstOrder)
		partialDependence = cbind(partialDependence, avg2ndOrder)		
		minDim = min(10, minDim2)
		varImportanceOverInteractions = vector()
		for (i in 1:minDim2)
		{
			idx = which(rownames(partialDependence)[i] == colnames(partialDependence))
			if (length(idx) > 0)  
			{	
				varImportanceOverInteractions[i] = 0.5*(avg1rstOrder[i] + avg2ndOrder[idx] + 
					2*partialDependence[i,idx])  
			}
			else 
			{	varImportanceOverInteractions[i] = avg1rstOrder[i]	}			
			names(varImportanceOverInteractions)[i] = rownames(partialDependence)[i]
		}				
		varImportanceOverInteractions = sort(varImportanceOverInteractions, decreasing = TRUE)		
		cat("\n\n2 - Local Variable importance")
		cat("\nVariables interactions (", minDim, " most important variables at first (columns) and second (rows) order) :", sep = "") 
		cat("\nFor each variable (at each order), its interaction with others is computed.\n\n")
		print(round(partialDependence[c(1:minDim, nrow(partialDependence)),],4))	
		cat("\n\nVariable Importance based on interactions (", minDim, " most important) :\n", sep = "")
		print(round(varImportanceOverInteractions[1:minDim]/sum(varImportanceOverInteractions), 4))
		if (!object$forest$regression)
		{
			cat("\nVariable importance over labels (", minDim, " most important variables conditionally to each label) :\n", sep = "")
			print(localVarImportance$classVariableImportance[1:minDim,])		
			cat("\n\nSee ...$localVariableImportance$obsVariableImportance to get variable importance for each observation.", sep = "")			
			cat("\n\nCall clusterAnalysis() function to get a more compact and complementary analysis.\n Type '?clusterAnalysis' for help.", sep = "")			
		}
		else
		{   cat("\n\nSee ...$localVariableImportance to get variable importance for each observation.")	}
		
		cat("\n\nCall partialDependenceOverResponses() function to get partial dependence over responses\nfor each variable. Type '?partialDependenceOverResponses' for help.\n")	
	}
	else
	{ stop("no variable importance defined in random uniform forest") }	
	importanceObject = list(globalVariableImportance = varImportance1, localVariableImportance = localVarImportance, 
		partialDependence = partialDependence, variableImportanceOverInteractions = varImportanceOverInteractions)		
	class(importanceObject) <- "importance"	
	importanceObject
}


partialImportance <- function(X, object, whichClass = NULL, threshold = NULL, thresholdDirection = c("low", "high"), border = NA, nLocalFeatures = 5)
{
	if (is.list(object$localVariableImportance))
	{
		Z = object$localVariableImportance$obsVariableImportance
		whichClassNames <- rm.string(names(object$localVariableImportance$class), "Class ")
		numericClassNames = as.numeric(as.factor(whichClassNames))		
		if (is.null(whichClass) ) { stop("Please provide a class.") }
		if (is.character(whichClass)) { whichClass = numericClassNames[which(whichClassNames == whichClass)] 	}
		whichClass2 = whichClassNames[which(numericClassNames == whichClass)]	
		idx = which(Z[,1] == whichClass)
		Z = Z[idx, ,drop = FALSE]		
		if (dim(Z)[1] <= 1) { stop("Not enough observations found using this class.") }
	}
	else
	{   
		Z = object$localVariableImportance	
		if (!is.null(threshold))
		{ 
			if (thresholdDirection == "low") { 	idx = which(Z[,1] <= threshold) }
			else {  idx = which(Z[,1] > threshold)  }
			Z = Z[idx,]
			
			if (dim(Z)[1] < 1) { stop("No observations found using this threshold.\n") }
		}
	}	
	Z = Z[,-1]
	idxLocalVar = grep("localVariable", colnames(X))
	idxPos = length(idxLocalVar)/2	
	countVars = sort(table(Z[,1:idxPos]), decreasing = TRUE)	
	obsObject = rmNA(countVars[1:nLocalFeatures]/sum(countVars))	
	XNames = colnames(X)
	par(las = 1)
	maxChar = floor(2 + max(nchar(XNames[as.numeric(names(sort(obsObject)))])))/2
	par(mar = c(5, maxChar + 2,4,2))
	barplot(sort(obsObject*100), horiz = TRUE, col = sort(heat.colors(length(obsObject)), decreasing = TRUE), 
	border = border, names.arg = XNames[as.numeric(names(sort(obsObject)))], xlab = "Relative influence (%)", 
		main = if (is.list(object$localVariableImportance)) { paste("Partial importance based on observations over class ", whichClass2, sep ="") }
		else 
		{
			if (!is.null(threshold)) 
			{	
				if (thresholdDirection == "low") { paste("Partial importance based on observations (with Response < ", round(threshold, 4), ")", sep ="") }
				else  { paste("Partial importance based on observations (with Response > ", round(threshold, 4), ")", sep ="") }					
			}
			else {	"Partial importance based on observations" }
		}
	)	
	cat("Relative influence: ", round(rmNA(sum(obsObject))*100, 2), "%\n", sep="")
	return(obsObject)
}

partialDependenceOverResponses <- function(Xtest, importanceObject, 
whichFeature = NULL, 
whichOrder = c("first", "second", "all"), 
outliersFilter = FALSE, 
plotting = TRUE, 
followIdx = FALSE, 
maxClasses = if (is.null(whichFeature)) { 10 } else { max(10, which.is.factor(Xtest[, whichFeature, drop = FALSE], count = TRUE)) }, 
bg = "lightgrey")
{
	FeatureValue = Response = Class = Observations = NULL
	if (!is.null(whichFeature))
	{	
		if (is.character(whichFeature)) { whichFeature = which(colnames(Xtest) == whichFeature) }		
		if (length(whichFeature) > 1)
		{ 
			whichFeature = whichFeature[1]
			cat("Only one variable can be computed at the same time\n")
		}
	}	
	if (whichOrder[1] == "first") 	{ idxOrder = 2	}	
	if (whichOrder[1] == "second") 	{ idxOrder = 3	}	
	if (whichOrder[1] == "all") 	
	{ 
		if (is.matrix(importanceObject$localVariableImportance))
		{	idxOrder = 2:length(grep("localVariableFrequency", colnames(importanceObject$localVariableImportance)))	}
		else
		{	idxOrder = 2:length(grep("localVariableFrequency", colnames(importanceObject$localVariableImportance$obsVariableImportance)))	}
	}	
	idx = list()
	if (is.matrix(importanceObject$localVariableImportance)) {	importanceObjectMatrix = importanceObject$localVariableImportance }
	else { importanceObjectMatrix = importanceObject$localVariableImportance$obsVariableImportance }	
	if (is.null(whichFeature))
	{	whichFeature = as.numeric(names(which.max(table(importanceObjectMatrix[,idxOrder[1]]))))	}			
	idx[[1]] = which(importanceObjectMatrix[,idxOrder[1]] == whichFeature)			
	if (length(idxOrder) > 1)
	{
		for (i in 1:length(idxOrder))
		{	idx[[i+1]] = which(importanceObjectMatrix[,1+idxOrder[i]]== whichFeature)	}
	}	
	partialDependenceMatrix = cbind(Xtest[unlist(idx), whichFeature], importanceObjectMatrix[unlist(idx),1], unlist(idx))
	partialDependenceMatrix = sortMatrix(partialDependenceMatrix, 1)	
	NAIdx = which(is.na(partialDependenceMatrix))	
	if (length(NAIdx) > 0) {  partialDependenceMatrix = partialDependenceMatrix[-NAIdx,] }	
	if (outliersFilter & (!is.factor(Xtest[,whichFeature])))
	{
		highOutlierIdx =  which(partialDependenceMatrix[,1] > quantile(partialDependenceMatrix[,1],0.95))
		lowOutlierIdx =  which(partialDependenceMatrix[,1] < quantile(partialDependenceMatrix[,1],0.05))
		if (length(highOutlierIdx) > 0 | length(lowOutlierIdx) > 0) 
		{	partialDependenceMatrix = partialDependenceMatrix[-c(lowOutlierIdx,highOutlierIdx),]	}
	}	
	if (is.vector(partialDependenceMatrix) )
	{ stop ("Not enough points to plot partial dependencies. Please increase order of interaction when computing importance.") }	
	if (dim(partialDependenceMatrix)[1] < 10)
	{ stop ("Not enough points to plot partial dependencies. Please increase order of interaction when computing importance.") }
	else
	{	
		idx = partialDependenceMatrix[,3]
	    partialDependenceMatrix = partialDependenceMatrix[,-3]
		flagFactor = 0
		smallLength = length(unique(Xtest[,whichFeature]))
		n = nrow(Xtest)		
		if ( ((smallLength < maxClasses) | is.factor(Xtest[,whichFeature])) & (smallLength/n <  1/5)) 
		{
			featureLevels = levels(as.factor(Xtest[,whichFeature]))
			testNumeric = is.numeric(as.numeric(featureLevels))
			if (testNumeric & length(rmNA(as.numeric(featureLevels))) > 0)
			{  flagFactor = 1 }
			else
			{	
				if (is.matrix(importanceObject$localVariableImportance))
				{  
					classFeature = unique(partialDependenceMatrix[,1])
					B = round(as.numeric(partialDependenceMatrix[,2]),4)		 
					A = as.numeric(factor2matrix(partialDependenceMatrix[,1, drop= FALSE]))
					partialDependenceMatrix = cbind(A,B)
					colnames(partialDependenceMatrix) = c("Class", "Response")
					flagFactor = 1 
					valueFeature = unique(A)
					
					referenceTab = cbind(classFeature, valueFeature)
					colnames(referenceTab) = c("category", "numeric value")
									
					cat("categorical values have been converted to numeric values :\n")
					print(referenceTab)
					cat("\n")
				}
				else
				{	partialDependenceMatrix[,1] = featureLevels[ as.numeric(as.factor(partialDependenceMatrix[,1]))] }
			}
		}
    }		
	if (plotting)
	{
		if (dim(partialDependenceMatrix)[1] < 1)
		{ stop ("Not enough points to plot partial dependencies. Please increase order of interaction when computing importance.") }	
		if (is.matrix(importanceObject$localVariableImportance))
		{
			## #require(ggplot2) || install.packages("ggplot2")
			if ( ((smallLength < maxClasses) | is.factor(Xtest[,whichFeature])) & (smallLength/n <  1/5)) 
		    {
				A = if (flagFactor) { as.factor(partialDependenceMatrix[,1]) } else { partialDependenceMatrix[,1] } 
				B = round(as.numeric(partialDependenceMatrix[,2]),4)				
				partialDependenceMatrix = data.frame(A , B)
				colnames(partialDependenceMatrix) = c("Class", "Response")
								
				plot(qplot(Class, Response, data = partialDependenceMatrix, geom = c("boxplot", "jitter"), 
				outlier.colour = "green", outlier.size = 2.5, fill= Class, main = "Partial dependence over predictor",
				xlab = colnames(Xtest)[whichFeature], ylab = "Response"))
			}
			else
			{			
				colnames(partialDependenceMatrix) = c("FeatureValue", "Response")
				partialDependenceMatrix = data.frame(partialDependenceMatrix)
			
				tt <- ggplot(partialDependenceMatrix, aes(x = FeatureValue, y = Response))
				plot(tt +  geom_point(colour = "lightblue") + stat_smooth(fill = "green", colour = "darkgreen", 
				size = 1) + 
				labs(title = "Partial dependence over predictor", x = colnames(Xtest)[whichFeature], y = "Response"))
			}
		}
		else
		{
			## #require(ggplot2) || install.packages("ggplot2")
			colnames(partialDependenceMatrix) = c("Observations", "Class")
			partialDependenceMatrix = data.frame(partialDependenceMatrix)
			variablesNames = unique(partialDependenceMatrix$Class)
			partialDependenceMatrix$Class = factor(partialDependenceMatrix$Class)
			levels(partialDependenceMatrix$Class) = colnames(importanceObject$localVariableImportance$classVariableImportance)[sort(variablesNames)]			
			if ( ((smallLength < maxClasses) | is.factor(Xtest[,whichFeature])) & (smallLength/n <  1/5))
			{   
				par(las=1)
				if (bg != "none") par(bg = bg)
				mosaicplot(t(table(partialDependenceMatrix)), color = sort(heat.colors(length(featureLevels)), 
			    decreasing = FALSE), border = NA, ylab = colnames(Xtest)[whichFeature],  xlab = "Class",
				main = "Partial dependence over predictor")
			}
			else
			{
				plot(qplot(Class, Observations, data = partialDependenceMatrix, geom = c("boxplot", "jitter"), 
				outlier.colour = "green", outlier.size = 2.5, fill= Class, main = "Partial dependence over predictor",
				xlab = "", ylab = colnames(Xtest)[whichFeature]))
			}
		}
	}
	else
	{
		if (!is.matrix(importanceObject$localVariableImportance))
		{	
			colnames(partialDependenceMatrix) = c("Observations", "Class")
			partialDependenceMatrix = data.frame(partialDependenceMatrix)
			variablesNames = unique(partialDependenceMatrix$Class)
			partialDependenceMatrix$Class = factor(partialDependenceMatrix$Class)
			levels(partialDependenceMatrix$Class) = colnames(importanceObject$localVariableImportance$classVariableImportance)[sort(variablesNames)]
		}
	}	
	if (followIdx)
	{	return(list(partialDependenceMatrix = partialDependenceMatrix, idx = as.numeric(idx) )) }
	else
	{  return(partialDependenceMatrix)  }
}	

partialDependenceBetweenPredictors <- function(Xtest, importanceObject, features, 
whichOrder = c("first", "second", "all"), 
perspective = FALSE, 
outliersFilter = FALSE, 
maxClasses = max(10, which.is.factor(Xtest[,features, drop = FALSE], count = TRUE)),
bg = "grey")
{
	Variable1 = Variable2 = SameClass = ..level.. = Response = NULL
	if (length(features) != 2) { stop("Please provide two features.") }	
	## #require(ggplot2)
	graphics.off()	
	pD1 <- partialDependenceOverResponses(Xtest, importanceObject, whichFeature = features[1], whichOrder = whichOrder, 
		outliersFilter = outliersFilter, plotting = FALSE, followIdx = TRUE, maxClasses = maxClasses)	
	pD2 <- partialDependenceOverResponses(Xtest, importanceObject, whichFeature = features[2], whichOrder = whichOrder, 
		outliersFilter = outliersFilter, plotting = FALSE, followIdx = TRUE, maxClasses = maxClasses)	
	sameIdx2 = find.idx(pD1$idx, pD2$idx)
	sameIdx1 = find.idx(pD2$idx, pD1$idx)
	minDim = length(sameIdx1)	
	if ( (minDim < 10) | (length(sameIdx2) < 10)) { stop("Not enough points. Please use option whichOrder = 'all'") }	
	pD1 = pD1$partialDependenceMatrix; pD2 = pD2$partialDependenceMatrix;
	pD11 = factor2matrix(pD1)[sameIdx1,]; pD22 = factor2matrix(pD2)[sameIdx2,]		
	if (!is.matrix(importanceObject$localVariableImportance))
	{
		minN = min(nrow(pD11), nrow(pD22))
		pD11 = pD11[1:minN,]
		pD22 = pD22[1:minN,]		
		idx = ifelse(pD11[,2] == pD22[,2], 1, 0)
		Xi  = cbind(pD11[which(idx == 1), 1], pD22[which(idx == 1), 1])
		Xj  = cbind(pD11[which(idx == 0), 1], pD22[which(idx == 0), 1])		
		if (!is.character(features[1]))
		{ 
			fName1 = which(colnames(Xtest) == colnames(Xtest)[features[1]])
			fName2 = which(colnames(Xtest) == colnames(Xtest)[features[2]])
			features = colnames(Xtest)[c(fName1, fName2)]
		}		
		Xi = as.data.frame(Xi); colnames(Xi) = c("Variable1", "Variable2")
		Xj = as.data.frame(Xj); colnames(Xj) = c("Variable1", "Variable2")		
		Xi = cbind(Xi, rep(1, nrow(Xi)))
		Xj = cbind(Xj, rep(0, nrow(Xj)))		
		colnames(Xj)[3] = colnames(Xi)[3] = "SameClass"
		X = rbind(Xi, Xj)
		X[,3] = ifelse(X[,3] == 1, TRUE, FALSE)
		smallLength = length(unique(Xtest[,features[1]]))
		n = nrow(Xtest)	
		Xnumeric = X
		ggplotFlag1 = ggplotFlag2 = 1
		options(warn = -1)
		if ( ((smallLength < maxClasses) | is.factor(Xtest[,features[1]])) & (smallLength/n <  1/5)) 
		{
			featureLevels = levels(factor(Xtest[,features[1]]))
			testNumeric = is.numeric(as.numeric(featureLevels))
			if (testNumeric & length(rmNA(as.numeric(featureLevels))) > 0)	{ 	ggplotFlag1 = 0 	}
			else  { X[,1] = featureLevels[X[,1]] }			
		}		
		smallLength = length(unique(Xtest[,features[2]]))
		if ( ((smallLength < maxClasses) | is.factor(Xtest[,features[2]])) & (smallLength/n < 1/5)) 
		{
			featureLevels = levels(factor(Xtest[,features[2]]))
			testNumeric = is.numeric(as.numeric(featureLevels))
			if (testNumeric & length(rmNA(as.numeric(featureLevels))) > 0)	{ 	ggplotFlag2 = 0 	}
			else  { X[,2] = featureLevels[X[,2]] }		
		}
		options(warn = 0)
		dev.new() 
		tt <- ggplot(X, aes(x = Variable1, y = Variable2, colour = SameClass))		
		if (ggplotFlag1 == 1)
		{
			plot(tt + geom_point(size = 2) 
				+ labs(title = "Dependence between predictors", x = features[1], y = features[2])	
				+ scale_colour_manual("Same class", values = c("red", "green") )
				+ theme(axis.text.x = element_text(angle = 60, hjust = 1))
			)
		}
		else
		{
			plot(tt + geom_point(size = 2) 
				+ labs(title = "Dependence between predictors", x = features[1], y = features[2])	
				+ scale_colour_manual("Same class", values = c("red", "green") )
			)
		}
		dev.new()
		cde1 <- geom_histogram(position = "fill", binwidth = diff(range(Xnumeric[,1]))/4, alpha = 7/10)
		cde2 <- geom_histogram(position = "fill", binwidth = diff(range(Xnumeric[,2]))/4, alpha = 7/10)
		tt1 <- ggplot(X, aes(x = Variable1, fill = SameClass))		
		if (ggplotFlag1 == 1)
		{
			plot(tt1 + cde1 
				+ labs(title = "Class distribution", x = features[1], y = "Frequency")
				+ scale_fill_manual(paste("Same class as", features[2]),values = c("red", "lightgreen"))
				+ theme(axis.text.x = element_text(angle = 60, hjust = 1))
			)
		}
		else
		{
			plot(tt1 + cde1 
				+ labs(title = "Class distribution", x = features[1], y = "Frequency")
				+ scale_fill_manual(paste("Same class as", features[2]),values = c("red", "lightgreen"))
			)
		}		
		dev.new()
		tt1 <- ggplot(X, aes(x = Variable2, fill = SameClass))		
		if (ggplotFlag2 == 1)
		{
			plot(tt1 + cde2 
				+ labs(title = "Class distribution", x = features[2], y = "Frequency")
				+ scale_fill_manual(paste("Same class as", features[1]),values = c("red", "lightgreen"))
				+ theme(axis.text.x = element_text(angle = 60, hjust = 1))
			)
		}
		else
		{
			plot(tt1 + cde2 
				+ labs(title = "Class distribution", x = features[2], y = "Frequency")
				+ scale_fill_manual(paste("Same class as", features[1]),values = c("red", "lightgreen"))
			)
		}		
		if ( (length(unique(X[,1])) == 1) | (length(unique(X[,2])) == 1) )
		{  cat("\nOne of the variable does not have variance. Heatmap is not plotted.\n") }
		else
		{
			dev.new()
			tt2 <- ggplot(Xnumeric, aes( x = Variable1, y = Variable2, z = SameClass))	
				try(plot(tt2 + stat_density2d(aes(fill = ..level.., alpha =..level..), geom = "polygon") 
				+ scale_fill_gradient2(low = "lightyellow", mid = "yellow", high = "red")
				+ labs(title = "Heatmap of dependence between predictors", x = features[1], y = features[2])
			), silent = TRUE)
		}
		colnames(X) = c(features, "Same class")
	}
	else
	{
		intervals = cut(c(pD11[,2], pD22[,2]), minDim, labels = FALSE)
		pD11 = cbind(pD11,intervals[1:minDim])		
		if (nrow(pD22) != length(rmNA(intervals[(minDim + 1):(2*minDim)])))
		{
			sameN = min(nrow(pD22),length(rmNA(intervals[(minDim + 1):(2*minDim)])))
			pD22 = cbind(pD22[1:sameN,], rmNA(intervals[(minDim + 1):(2*minDim)])[1:sameN])
		}
		else
		{ 	pD22 = cbind(pD22, rmNA(intervals[(minDim + 1):(2*minDim)]))	}		
		minN = min(nrow(pD11), nrow(pD22))		
		Xi = sortMatrix(pD11,3)[1:minN,]
		Xj = sortMatrix(pD22,3)[1:minN,]
		Z  = (Xi[,2] + Xj[,2])/2		
		if (!is.character(features[1]))
		{ 
			fName1 = which(colnames(Xtest) == colnames(Xtest)[features[1]])
			fName2 = which(colnames(Xtest) == colnames(Xtest)[features[2]])
			features = colnames(Xtest)[c(fName1, fName2)]
		}		
		dev.new()
		XiXj = as.data.frame(cbind(Xi[,1], Xj[,1], Z)); colnames(XiXj) = c("Variable1", "Variable2", "Response")		
		tt <- ggplot(XiXj, aes(x = Variable1, y = Variable2, z = Response))
		ttMore <- tt + stat_density2d(aes(fill = ..level.., alpha =..level..), geom = "polygon") +
			scale_fill_gradient2(low = "lightyellow", mid = "yellow", high = "red") +
			labs(title = "Local Heatmap of dependence (with frequency and intensity of Response)", 
			x = features[1], y = features[2])
		try(plot(ttMore), silent= TRUE)		
		dev.new()
		fourQuantilesCut = cut(Z,4)
		XiXj[,3] = fourQuantilesCut		
		dataCuts = table(fourQuantilesCut)		
		if (length(which(dataCuts < 5)) > 0)
		{
			lowDataCuts = names(dataCuts[which(dataCuts < 5)])
			rmIdx = find.idx(lowDataCuts, XiXj[,3])			
			XiXj = XiXj[,-3]
			Z = Z[-rmIdx]
			fourQuantilesCut = cut(Z, 4)
			XiXj = cbind(XiXj[-rmIdx,], fourQuantilesCut)		
			dataCuts = table(fourQuantilesCut)
			colnames(XiXj)[3] = "Response"
		}
		X = cbind(XiXj, Z)
		colnames(X) = c(features, "Response in four quantile intervals", "Response")		
	    tt1 <- ggplot(XiXj, aes(x = Variable1, y = Variable2, z = Response)) 
		try(plot(tt1 + stat_density2d(aes(fill = ..level.., alpha =..level..), geom = "polygon") 
			+  scale_fill_gradient2(low = "lightyellow", mid = "yellow", high = "red")	
			+ labs(title = "Global Heatmap of dependence (with intensity of Response)", x = features[1], 
			y = features[2]) 
		), silent = TRUE)				
		dev.new()
		try(plot(tt + geom_point(aes(colour = Response, size = Response))
		    +  stat_smooth(fill = "lightgrey", colour = "grey", size = 1) 
			+ labs(title = "Dependence between predictors", x = features[1], y = features[2])		
			+ scale_colour_gradient2(low = "blue", mid = "green", high = "red") 	
		), silent = TRUE)
	}	
	#if (.Platform$OS.type == "windows") { arrangeWindows("vertical", preserve = FALSE) }
	idxf1 = which(colnames(importanceObject$partialDependence) ==  features[1])
	idxf2 = which(rownames(importanceObject$partialDependence) ==  features[2])	
	if ( (length(idxf1) > 0) &  (length(idxf2) > 0) ) 
	{
		np = dim(importanceObject$partialDependence)
		cat("\nLevel of interactions between ", features[1], " and " , features[2], " at first order: ", 
		round(importanceObject$partialDependence[idxf2, idxf1],4), "\n", "(", round(round(importanceObject$partialDependence[idxf2, idxf1],4)/max(importanceObject$partialDependence[-np[1], -np[2]])*100,2), "% of the feature(s) with maximum level", ")\n",  sep="")
	}
	else
	{  cat("\nFeatures do not appear to have strong co-influence on the response.\n") }	
	idxf1 = which(rownames(importanceObject$partialDependence) == features[1])	
	if (length(idxf1) > 0)	
	{ 
		idxf2 = which(colnames(importanceObject$partialDependence) ==  features[2]) 
		if (length(idxf2) > 0)	
		{	
			cat("Level of interactions between ", features[1], " and ", features[2], " at second order: ", 
			round(importanceObject$partialDependence[idxf1, idxf2],4), "\n", "(", round(round(importanceObject$partialDependence[idxf1, idxf2],4)/max(importanceObject$partialDependence[-np[1], -np[2]])*100,2), "% of the feature(s) with maximum level", ")\n",  sep="")
		}
	}
	
	if(!is.matrix(importanceObject$localVariableImportance))
	{	
		cat("\nClass distribution : for a variable of the pair, displays the estimated probability\nthat the considered variable has the same class than the other. If same class tends to be TRUE\nthen the variable has possibly an influence on the other (for the considered category or values) when predicting a label.\n\n")
		cat("Dependence : for the pair of variables, displays the shape of their dependence\nand the estimated agreement in predicting the same class, for the values that define dependence.\nIn case of categorical variables, cross-tabulation is used.\n\n")
		cat("Heatmap : for the pair of variables, displays the area where the dependence is the most effective.\nThe darker the colour, the stronger is the dependence.\n")
		cat("\nFrom the pair of variables, the one that dominates is, possibly, the one\nthat is the most discriminant one (looking 'Global variable Importance') and/or the one\nthat has the higher level of interactions(looking 'Variable Importance based on interactions').\n")
	}
	else
	{
		cat("\nDependence : for the pair of variables, displays the shape of their dependence\nand the predicted value (on average) of the response for the values taken by the pair.\n\n")
		cat("Heatmap : for the pair of variables, displays the area where the dependence is the most effective.\nThe darker the colour, the stronger is the dependence. First plot focuses on the intensity of the response\n, while the second considers both frequency (number of close predictions) and intensity.\n\n")	
	}		
	cat("\nPlease use the R menu to tile vertically windows in order to see all plots.\n\n")
	if (perspective & is.matrix(importanceObject$localVariableImportance))
	{
		dev.new()
		gridSize = 100; gridLag = 100
		x = XiXj[,1]
		y = XiXj[,2]
		z = Z			
		n = length(x)		
		xyz = cbind(x,y,z)
		if (n > 5*gridSize)
		{
			sampleIdx = sample(n, 500)
			xyz = xyz[sampleIdx,]
			n = length(sampleIdx)
		}
		xyz.byX = sortMatrix(xyz,1)
		xyz.byY = sortMatrix(xyz,2)			
		newX = xyz.byX[,1]
		newY = xyz.byY[,2]		
		dummyForRepX = dummyForRepY = rep(0,n)
		for (i in 2:n)
		{ 
			if (newX[i] == newX[i-1]) {	dummyForRepX[i] = 1 }
			if (newY[i] == newY[i-1]) {	dummyForRepY[i] = 1 }
		}			
		newIdx = which(dummyForRepY == 0 & dummyForRepX == 0)
		newX = newX[newIdx]
		newY = newY[newIdx]			
		if ( (gridSize + gridLag) > length(newX))
		{
			interp.1 = seq(min(newX) + 0.1, max(newX) - 0.1,length = gridSize + gridLag - length(newX))
			interp.2 = seq(min(newY) + 0.1, max(newY) - 0.1,length = gridSize + gridLag- length(newY))				
			newX = sort(c(newX, interp.1))
			newY = sort(c(newY, interp.2))
		}			
		duplicatesX = duplicated(newX)
		duplicatesY = duplicated(newY)			
		if (sum(duplicatesX) > 0)
		{
			newX = newX[!duplicatesX]
			newY = newY[!duplicatesX]
		}			
		if (sum(duplicatesY) > 0)
		{
			newY = newY[!duplicatesY]
			newX = newX[!duplicatesY]
		}			
		newXYZ = cbind(newX, newY, rep(NA, length(newX)))
		proxyM = fillNA2.randomUniformForest(rbind(xyz, newXYZ), nodesize = 2)						
		xyz.dim =  dim(xyz)
		nn = length(newX)						
		newZ = matrix(NA, nn, nn)			
		for (i in 1:nn)
		{
			for (j in 1:nn)
			{	newZ[i,j] = mean(proxyM[which(newX[i] == proxyM[,1] | newY[j] == proxyM[,2]), 3])	}
		}											
		L.smoothNewZ = t(apply(newZ, 1, function(Z) lagFunction(Z, lag = gridLag, FUN  = mean, inRange = TRUE)))
		C.smoothNewZ = apply(newZ, 2, function(Z) lagFunction(Z, lag = gridLag, FUN  = mean, inRange = TRUE))			
		smoothIdx = round(seq(1, nrow(L.smoothNewZ), length = gridLag),0)
		newX = newX[-smoothIdx]
		newY = newY[-smoothIdx]
		C.smoothNewZ = C.smoothNewZ[,-smoothIdx]
		L.smoothNewZ = L.smoothNewZ[-smoothIdx,]
		newZ = 0.5*(C.smoothNewZ + L.smoothNewZ)		
		highOutlierIdx2 <- apply(newZ, 2, function(Z) which(Z >= quantile(Z, 0.975)))
		lowOutlierIdx2 <- apply(newZ, 2, function(Z) which(Z <= quantile(Z, 0.025)))
		ouliersIdx2 <- c(lowOutlierIdx2, highOutlierIdx2)
		if (length(ouliersIdx2) > 0) 
		{	
			newZ = newZ[-ouliersIdx2, -ouliersIdx2]  	
			newX = newX[-ouliersIdx2]
			newY = newY[-ouliersIdx2]
		}
		rm(lowOutlierIdx2);
		rm(highOutlierIdx2);
		highOutlierIdx2 <- apply(newZ, 1, function(Z) which(Z >= quantile(Z, 0.975)))
		lowOutlierIdx2 <-  apply(newZ, 1, function(Z) which(Z <= quantile(Z, 0.025)))
		ouliersIdx2 <- c(lowOutlierIdx2, highOutlierIdx2)
		if (length(ouliersIdx2) > 0) 
		{	
			newZ = newZ[-ouliersIdx2, -ouliersIdx2]  	
			newX = newX[-ouliersIdx2]
			newY = newY[-ouliersIdx2]
		}			
		nNewZ = nrow(newZ)
		if (nNewZ > gridSize)
		{
			sampleIdx = sort(sample(nNewZ, gridSize))
			newX = newX[sampleIdx]
			newY = newY[sampleIdx]
			newZ = newZ[sampleIdx, sampleIdx]
		}
		if (bg != "none") par(bg = bg)
		flag = endCondition = 0
		lastANSWER = -40
		newX = as.matrix(newX)
		newY = as.matrix(newY)
		newZ = as.matrix(newZ)
		while (!endCondition)
		{
			if (!flag)
			{
				try(perspWithcol(newX, newY, newZ, heat.colors, nrow(newZ), theta = -40, phi = 20, xlab = features[1], ylab = features[2], zlab = "Response", main = "Dependence between predictors and effect over Response", ticktype = "detailed", box = TRUE, expand = 0.5, shade = 0.15), silent = FALSE)
			}			
			ANSWER <- readline(cat("To get another view of 3D representation\nplease give a number between -180 and 180 (default one = -40).\nType 'b' to remove border.\nTo see animation please type 'a'.\nType Escape to leave :\n"))
			if ( is.numeric(as.numeric(ANSWER)) & !is.na(as.numeric(ANSWER)) )
			{  
				flag = 1
				try(perspWithcol(newX, newY, newZ, heat.colors, nrow(newZ), theta = ANSWER, phi = 20, xlab = features[1], ylab = features[2], zlab = "Response", main = "Dependence between predictors and effect over Response", ticktype = "detailed", box = TRUE, expand = 0.5, shade = 0.15), 
				silent = FALSE)
				lastANSWER = ANSWER
			}
			else
			{ 
				if (as.character(ANSWER) == "b")
				{
					flag = 1
					try(perspWithcol(newX, newY, newZ, heat.colors, nrow(newZ), theta = lastANSWER, phi = 20, xlab = features[1], ylab = features[2], zlab = "Response", main = "Dependence between predictors and effect over Response", ticktype = "detailed", box = TRUE, expand = 0.5, border = NA, shade = 0.15), silent = FALSE)
				}
				else
				{
					if ( as.character(ANSWER) == "a")
					{
						flag = 1
						nn = nrow(newZ)
						circle = -180:180
						for (i in circle)
						{
							try(perspWithcol(newX, newY, newZ, heat.colors, nn, theta = i, phi = 20, 
							xlab = features[1], ylab = features[2], zlab = "Response", main = "Dependence between predictors and effect over Response", ticktype = "detailed", box = TRUE, border = NA, expand = 0.5, shade = 0.15), silent = FALSE)	
						}
					}				
					else				
					{ endCondition = 1	}
				}
			}
		}			
	}		
	return(X)
}

twoColumnsImportance <- function(importanceObjectMatrix)
{
	idx = length(grep("localVariableFrequency", colnames(importanceObjectMatrix))) - 1
	tmpImportanceObjectMatrix = importanceObjectMatrix[,1:2]
	if (idx > 0)
	{
		for (i in 1:idx)
		{	tmpImportanceObjectMatrix = rbind(tmpImportanceObjectMatrix, importanceObjectMatrix[,c(1,2+idx[i])])	}
	}
	return(tmpImportanceObjectMatrix)
}

plot.importance <- function(x, 
nGlobalFeatures = 30, 
nLocalFeatures = 5, 
Xtest = NULL, 
whichFeature = NULL, 
whichOrder = if (ncol(x$globalVariableImportance) > 5) 
	{  if (nrow(x$localVariableImportance$obsVariableImportance) > 1000) "first" else "all" } 
	else
	{  if (nrow(x$localVariableImportance) > 1000) "first" else "all" }, 
outliersFilter = FALSE, 
formulaInput = NULL, 
border = NA, 
...)
{
	cat("\nPlease use the R menu to tile vertically windows in order to see all plots.\n")
	object <- x
	if (nrow(x$partialDependence) < 3)
	{
		stop("Not enough (or no) interactions between variables. Please use partialDependenceOverResponses( )\n and partialDependenceBetweenPredictors( ) functions to deeper assess importance.") 
	}
	Variable = Response = NULL
	if (!is.null(Xtest))
	{
		if (!is.null(formulaInput))
		{
			mf <- model.frame(formula = formulaInput, data = as.data.frame(Xtest))
			Xtest <- model.matrix(attr(mf, "terms"), data = mf)[,-1]
			cat("Note: please note that categorical variables have lost their original values when using formula.\nIt is strongly recommended to not use formula if one wants to assess importance\n \n.")
		}
		else
		{		
			matchNA = (length(which(is.na(Xtest))) > 0)				
			if (matchNA) 
			{ 
				cat("NA found in data. Fast imputation (means) is used for missing values\n")
				Xtest <- na.impute(Xtest)
			}
		}
	}		
	maxVar = nGlobalFeatures
	maxVar2 = nLocalFeatures
	graphics.off()
	varImportance = object$globalVariableImportance
	n = nrow(varImportance)	
	if (n > maxVar) { varImportance = varImportance[1:maxVar,] }
	else { maxVar = n}	
	par(las = 1)
	maxChar = if (nLocalFeatures >= 10) { 11 } 
	else { floor(2 + max(nchar(as.character(object$globalVariableImportance[,1])))/2) }
	par(mar = c(5, maxChar + 1,4,2))
	barplot(varImportance[maxVar:1,"percent.importance"], horiz = TRUE, col = sort(heat.colors(maxVar), decreasing = TRUE), border = border,
	names.arg = varImportance[maxVar:1,"variables"], xlab = "Relative influence (%)", main = "Variable importance based on information gain")
	abline(v = 100/n, col = 'grey')	
	dev.new()
	par(las=1)
    par(mar = c(5,4,4,2))
	nbFeatures = ncol(object$partialDependence)
	newNbFeatures = min(maxVar2, nbFeatures -1)	
	if (newNbFeatures < (nbFeatures - 1))
	{
		OthersVariablesCol = colSums(object$partialDependence[(newNbFeatures+1):(nbFeatures -1), -nbFeatures, drop = FALSE])[1:newNbFeatures]
		OthersVariablesRow = rowSums(object$partialDependence[-nbFeatures, (newNbFeatures+1):(nbFeatures -1), drop = FALSE])[1:newNbFeatures]
		corner = mean(c(OthersVariablesCol, OthersVariablesRow))
		newPartialDependence = object$partialDependence[1:newNbFeatures,1:newNbFeatures]
		newPartialDependence =  rbind(cbind(newPartialDependence, OthersVariablesRow), c(OthersVariablesCol, corner))
		colnames(newPartialDependence)[ncol(newPartialDependence)] = "Other features"
		rownames(newPartialDependence)[nrow(newPartialDependence)] = "Other features"
		mosaicplot(t(newPartialDependence), color = sort(heat.colors(newNbFeatures + 1), decreasing = FALSE), 
		main = "Variables interactions over observations", ylab = "Most important variables at 2nd order", 
		xlab = "Most important variables at 1rst order", las = ifelse(maxChar > 10, 2,1), border = border)
	}
	else
	{
		mosaicplot(t(object$partialDependence[1:newNbFeatures,1:newNbFeatures]), color = sort(heat.colors(newNbFeatures), decreasing = FALSE), 
		las = ifelse(maxChar > 10, 2, 1), main = "Variables interactions over observations", ylab = "Most important variables at 2nd order", xlab = "Most important variables at 1rst order", border = border)
	}	
	dev.new()
	par(las=1)
	par(mar=c(5,maxChar + 1,4,2))
	nbFeatures2 = min(maxVar2, length(object$variableImportanceOverInteractions))
	importanceOverInteractions = sort(object$variableImportanceOverInteractions, decreasing = TRUE)/sum(object$variableImportanceOverInteractions)*100
	barplot(importanceOverInteractions[nbFeatures2:1], horiz = TRUE, col = sort(heat.colors(nbFeatures2), decreasing = TRUE), 
		names.arg = names(importanceOverInteractions)[nbFeatures2:1], xlab = "Relative influence (%)", 
		main = "Variable importance based on interactions", border = border)
	abline(v = 100/n, col = 'grey')
	if (!is.matrix(object$localVariableImportance))
	{
		dev.new()
		par(las=1)
	    par(mar = c(5, maxChar + 1,4,2))
		nbFeatures3 = min(nbFeatures2, nrow(object$localVariableImportance$classVariableImportance))
		mosaicplot(t(object$localVariableImportance$classVariableImportance[1:nbFeatures3,,drop = FALSE]), 
			color = sort(heat.colors(nbFeatures3),
			decreasing = FALSE), main = "Variable importance over labels", border = border)
	}
	else
	{
		dev.new()
		## #require(ggplot2) || install.packages("ggplot2")
		ggData = twoColumnsImportance(object$localVariableImportance)
		mostImportantFeatures = as.numeric(names(sort(table(ggData[,2]), decreasing = TRUE)))	
		if (length(unique(ggData[,2])) > maxVar2)
		{	
			mostImportantFeatures = mostImportantFeatures[1:maxVar2]	
			ggData = ggData[find.idx(mostImportantFeatures, ggData[,2], sorting = FALSE),]
		}
		if (is.null(Xtest)) {	textX = paste("V", mostImportantFeatures, collapse = " ", sep="")	}
		else {	textX = paste( colnames(Xtest)[mostImportantFeatures], collapse = ", ", sep="")	}
		textX = paste("Index of most important variables [", textX, "]")
		colnames(ggData)[1] = "Response"
		colnames(ggData)[2] = "Variable"
		ggData = data.frame(ggData)
		if (!is.null(Xtest)) { 	ggData[,"Variable"] = colnames(Xtest)[ggData[,"Variable"]] 	}
		ggData[,"Variable"] = as.factor(ggData[,"Variable"])
		if (nrow(ggData) > 1000) 
		{ 
			randomSample = sample(nrow(ggData), 1000) 
			gg <- qplot(Variable, Response, data = ggData[randomSample,], geom = c("boxplot", "jitter"), 
			colour = Response, outlier.colour = "green", outlier.size = 1.5, fill = Variable)
			cat("1000 observations have been randomly sampled to plot 'Dependence on most important predictors'.\n")
		}
		else
		{	
			gg <- qplot(Variable, Response, data = ggData, geom = c("boxplot", "jitter"), colour = Response, outlier.colour = "green", outlier.size = 1.5, fill = Variable)	
		}		
		if (maxChar > 10)
		{
			plot(gg + labs(x ="", y = "Response", title = "Dependence on most important predictors") + 
			theme(axis.text.x = element_text(angle = 60, hjust = 1)))			
		}
		else
		{	
			plot(gg + labs(x ="", y = "Response", title = "Dependence on most important predictors"))
		}
	}	
	if (is.null(Xtest))	{	stop("partial dependence between response and predictor can not be computed without data") }
	else
	{
		endCondition = 0
		dev.new()
		pD <- partialDependenceOverResponses(Xtest, object, whichFeature = whichFeature, whichOrder = whichOrder, outliersFilter = outliersFilter,...)
		#if (.Platform$OS.type == "windows") { arrangeWindows("vertical", preserve = FALSE) }
		idxMostimportant = rmNA(match(names(object$variableImportanceOverInteractions), colnames(Xtest)))[1:nbFeatures2]
		mostimportantFeatures = colnames(Xtest)[idxMostimportant]
		while (!endCondition)
		{
			ANSWER <- readline(cat("To get partial dependence of most important features\ngive a column number\namong", idxMostimportant,
			"\n(", mostimportantFeatures, ")\nPress escape to quit\n"))
			if ( is.numeric(as.numeric(ANSWER)) & !is.na(as.numeric(ANSWER)) )
			{  
				whichFeature = as.numeric(ANSWER)
				if (whichFeature %in% idxMostimportant)
				{	
					pD <- partialDependenceOverResponses(Xtest, object, whichFeature = whichFeature, whichOrder = whichOrder, outliersFilter = outliersFilter,...)	
				}
				else
				{   stop("Please provide column index among most important. Partial Dependence can not be computed.") }				
			}
			else
			{	endCondition = 1	}
		}
	}	
}


print.importance <- function(x,...)
{
	object <- x
	minDim = min(10,length(object$variableImportanceOverInteractions))
	# global importance
	if (!is.matrix(object$localVariableImportance))
	{
		cat("\n1 - Global Variable importance (", minDim, " most important based on information gain) :\n", sep = "")
		cat("Note: most predictive features are ordered by 'score' and plotted. Most discriminant ones\nshould also be taken into account by looking 'class' and 'class.frequency'.\n\n")
	}
	else
	{  
		cat("\n1 - Global Variable importance (", minDim, " most important based on 'Lp' distance) :\n", 
		sep = "") 
	}
	print(object$globalVariableImportance[1:minDim,])	
	
	cat("\n\n2 - Local Variable importance")
	cat("\nVariables interactions (", minDim, " most important variables at first (columns) and second (rows) order) :", sep = "")
	cat("\nFor each variable (at each order), its interaction with others is computed.\n\n")
	#interactions	
	print(round(object$partialDependence[c(1:minDim, nrow(object$partialDependence)),], 2))	
	#local importance
	cat("\n\nVariable Importance based on interactions (", minDim, " most important) :\n", sep = "")
	print(round(object$variableImportanceOverInteractions/sum(object$variableImportanceOverInteractions),2)[1:minDim])	
	#local importance over labels
	if (!is.matrix(object$localVariableImportance))
	{
		cat("\nVariable importance over labels (", minDim, 
		" most important variables conditionally to each label) :\n", sep = "")
		print(object$localVariableImportance$classVariableImportance[1:minDim,])
	
		cat("\n\nSee ...$localVariableImportance$obsVariableImportance to get variable importance for each observation.", sep = "")
		
		cat("\n\nCall clusterAnalysis() function to get a more compact and complementary analysis.\n Type '?clusterAnalysis' for help.", sep = "")	
	}
	else 
	{   cat("\n\nSee ...$localVariableImportance to get variable importance for each observation.")	}	
	# partial dependence (response vs predictor)
	cat("\n\nCall partialDependenceOverResponses() function to get partial dependence over responses\nfor each variable. Type '?partialDependenceOverResponses' for help.\n")	
}

combineRUFObjects <- function(rUF1, rUF2)
{
	rUF1 <- filter.object(rUF1)
	rUF2 <- filter.object(rUF2)
	return(onlineCombineRUF(rUF1, rUF2))	
}

rankingTrainData <- function(trainData = NULL, trainLabels = NULL,  testData = NULL, testLabels = NULL, ntree = 100,  thresholdScore = 2/3, nTimes = 2, ...)
{	
	score = tmpScore = rep(0, nrow(testData))	
	classes = unique(trainLabels)
	majorityClass = modX(trainLabels)	
	i = 1;	rmIdx = NULL; idx = 1:nrow(testData)
	while  ( (i <= nTimes)  & (length(testData[,1]) >= 1) )
	{
		rUF <- randomUniformForestCore(trainData, trainLabels = as.factor(trainLabels), ntree = ntree, use.OOB = FALSE, rf.overSampling = -0.75, rf.targetClass = majorityClass, rf.treeSubsampleRate = 2/3)
		predictRUF <- randomUniformForestCore.predict(rUF, testData)	
		tmpScore = tmpScore + apply(predictRUF$all.votes, 1, function(Z)  length(which(Z != majorityClass)))
		score[idx] = score[idx] + tmpScore								
		rmIdx = which(tmpScore < (thresholdScore*ntree))		
		if (length(rmIdx) > 0)  {  testData = testData[-rmIdx,];  idx = idx[-rmIdx]; tmpScore = tmpScore[-rmIdx];   }
		i = i + 1
	}	
	return(score)
}
		
plotTreeCore <- function(treeStruct, rowNum = 1, height.increment = 1)
{
  if ( (treeStruct[rowNum, "status"] == -1) )
  {
    treeGraphStruct <- list()
    attr(treeGraphStruct, "members") <- 1
    attr(treeGraphStruct, "height") <- 0
    attr(treeGraphStruct, "label") <- if (treeStruct[rowNum,"prediction"] == 0) { "next node" } else 
		{ 	
			if (is.numeric(treeStruct[rowNum,"prediction"]))	
			{ round(treeStruct[rowNum,"prediction"],2) } else { treeStruct[rowNum,"prediction"]	}	
		}
	attr(treeGraphStruct, "leaf") <- TRUE
  }
  else
  {
	left <- plotTreeCore(treeStruct, treeStruct[rowNum, "left.daughter"], height.increment) #,  incDepth+1
    right <- plotTreeCore(treeStruct, treeStruct[rowNum, "right.daughter"], height.increment)#,  incDepth+1
    treeGraphStruct <- list(left,right)
    attr(treeGraphStruct, "members") <- attr(left, "members") + attr(right,"members")
    attr(treeGraphStruct,"height") <- max(attr(left, "height"),attr(right, "height")) + height.increment
    attr(treeGraphStruct, "leaf") <- FALSE	
	if (rowNum != 1)
	{     attr(treeGraphStruct, "edgetext") <- paste(treeStruct[rowNum, "split.var"] , " > " , round(treeStruct[rowNum, "split.point"],2), " ?", sep ="")  }
	else
	{	  
		attr(treeGraphStruct, "edgetext") <- paste(".[no].  .", 
		treeStruct[rowNum, "split.var"] , " > " , round(treeStruct[rowNum, "split.point"],2), ".                                                   .[yes]", sep="") 
	}
  }
  class(treeGraphStruct) <- "dendrogram"
  return(treeGraphStruct)
}

plotTreeCore2 <- function(treeStruct, rowNum = 1, height.increment = 1,  maxDepth = 100 )
{
  if ((treeStruct[rowNum, "status"] == -1) | (rowNum > maxDepth))
  {
    treeGraphStruct <- list()
    attr(treeGraphStruct, "members") <- 1
    attr(treeGraphStruct, "height") <- 0
    attr(treeGraphStruct, "label") <- if (treeStruct[rowNum,"status"] == 1) { "next node" } else 
		{ 	
			if (is.numeric(treeStruct[rowNum,"prediction"]))	
			{ round(treeStruct[rowNum,"prediction"],2) } else { treeStruct[rowNum,"prediction"]	}	
		}
	attr(treeGraphStruct, "leaf") <- TRUE
  }
  else
  {
    left <- plotTreeCore2(treeStruct, treeStruct[rowNum, "left.daughter"], height.increment) #,  incDepth+1
    right <- plotTreeCore2(treeStruct, treeStruct[rowNum, "right.daughter"], height.increment)#,  incDepth+1
    treeGraphStruct <- list(left,right)
    attr(treeGraphStruct, "members") <- attr(left, "members") + attr(right,"members")
    attr(treeGraphStruct,"height") <- max(attr(left, "height"),attr(right, "height")) + height.increment
    attr(treeGraphStruct, "leaf") <- FALSE	
	if (rowNum != 1)
	{     attr(treeGraphStruct, "edgetext") <- paste(treeStruct[rowNum, "split.var"] , " > " , round(treeStruct[rowNum, "split.point"],2), "?", sep ="")  }
	else
	{	  
		attr(treeGraphStruct, "edgetext") <- paste(treeStruct[rowNum, "split.var"] , " > " , round(treeStruct[rowNum, "split.point"],2), "?", ".                                     .", sep="") 
	}
  }
  class(treeGraphStruct) <- "dendrogram"
  return(treeGraphStruct)
}

plotTree <- function(treeStruct, rowNum = 1, height.increment = 1, maxDepth = 100, fullTree = FALSE, xlim = NULL, ylim= NULL, center = TRUE)
{
	if (fullTree)
	{	drawTreeStruct <- plotTreeCore(treeStruct, rowNum = rowNum , height.increment = height.increment)	}
	else
	{    drawTreeStruct <- plotTreeCore2(treeStruct, rowNum = rowNum , height.increment = height.increment,  maxDepth = maxDepth)	}	
	nP <- list(col = 3:2, cex = c(2.0, 0.75), pch = 21:22, bg =  c("light blue", "pink"), lab.cex = 0.75, lab.col = "tomato")	
	if (is.null(xlim))
	{	
		if (is.null(ylim)) { ylim = c(0,8) }
		plot(drawTreeStruct, center = center, leaflab ='perpendicular', edgePar = list(t.cex = 2/3, p.col = NA, p.lty = 0, lty =c( 2,5), 
			col = c("purple", "red"), lwd = 1.5), nodePar = nP, ylab = "Tree depth", xlab = "Predictions", 
			xlim = c(0, min(30,floor(nrow(treeStruct)/2))), ylim = ylim)
	}
	else
	{
		if (is.null(ylim)) { ylim = c(0,8) }
		plot(drawTreeStruct, center = center, leaflab ='perpendicular', edgePar = list(t.cex = 2/3, p.col = NA, p.lty = 0, lty =c( 2,5), 
			col = c("purple", "red"), lwd = 1.5), nodePar = nP, ylab = "Tree depth", xlab = "Predictions", xlim = xlim, ylim = ylim )
	}
}

fillNA2.randomUniformForest <- function(X, Y = NULL, ntree = 100, mtry = 1, nodesize = 10, 
categoricalvariablesidx = NULL, 
NAgrep = "", 
maxClasses = floor(0.01*min(3000, nrow(X))+2), 
threads = "auto",
...)
{	
    i = NULL 
	n <- nrow(X)
	X <- fillVariablesNames(X)	
	if (exists("categorical", inherits = FALSE)) { categoricalvariablesidx = categorical }
	else { categorical = NULL  }
	if (!is.null(Y)) { trueY = Y }
	trueX = X
	limitSize = if (n > 2000) { 50 } else { max(nodesize, 10) } 
	if (!is.matrix(X))
	{	
		flag = TRUE
		X.factors <- which.is.factor(X, maxClasses = maxClasses)
		X <- NAfactor2matrix(X, toGrep = NAgrep)
	}
	else
	{
		flag = FALSE
		X.factors = rep(0,ncol(X))
	}	
	NAIdx = which(is.na(X), arr.ind = TRUE)	
	if (dim(NAIdx)[1] == 0) { stop("No missing values in data. Please use NAgrep option to specify missing values.\nFunction checks for NA in data and for the string given by NAgrep. Please also check if string does not contain space character") }	
	processedFeatures = unique(NAIdx[,2])
	nbNA = table(NAIdx[,2])
	fullNAFeatures = which(nbNA == n)
	fullNAFeaturesLength = length(fullNAFeatures)
	if (fullNAFeaturesLength == length(processedFeatures))
	{	stop("All features have only NA in their values.\n") }
	else
	{	
		if (fullNAFeaturesLength > 0) 
		{	processedFeatures = processedFeatures[-fullNAFeatures] }
	}
	nFeatures = length(processedFeatures)	
	idx <- lapply(processedFeatures, function(Z) NAIdx[which(NAIdx[,2] == Z),1]) 
	validIdx <- sapply(idx, function(Z) (length(Z) >  limitSize))
	processedFeatures <- processedFeatures[which(validIdx == TRUE)]
	nFeatures = length(processedFeatures)
	invalidIdx <- which(validIdx == FALSE)
	if (length(invalidIdx) > 0)	{	idx <- rm.InAList(idx, invalidIdx)	}		
	if (nFeatures == 0)
	{	
		cat("Not enough missing values. Rough imputation is done. Please lower 'nodesize' value to increase accuracy.\n")
		return(na.replace(trueX, fast = TRUE))	
	}
	else
	{
		if (!is.null(Y)){  X = cbind(X,Y)	}
		X <- fillVariablesNames(X)
		X <- na.replace(X, fast = TRUE)			
		{
			## #require(parallel)	
			max_threads = min(detectCores(),4)		
			if (threads == "auto")
			{	
				if (max_threads == 2) { threads = min(max_threads, nFeatures) }
				else {	threads = max(1, max_threads)  }
			}
			else
			{
				if (max_threads < threads) 
				{	cat("Note: number of threads is higher than number of logical threads in this computer.\n") }
			}			
			threads = min(nFeatures, max_threads)
			## #require(doParallel)			
			Cl <- makePSOCKcluster(threads, type = "SOCK")
			registerDoParallel(Cl)		
			chunkSize <- ceiling(nFeatures/getDoParWorkers())
			smpopts  <- list(chunkSize = chunkSize)
		}				
		export = c("randomUniformForest.default", "rUniformForest.big", "randomUniformForestCore.big", "randomUniformForestCore", "predict.randomUniformForest", "rUniformForestPredict", "uniformDecisionTree", "CheckSameValuesInAllAttributes", "CheckSameValuesInLabels", "fullNode", "genericNode", "leafNode", "filter.object", "filter.forest", "randomUniformForestCore.predict", "onlineClassify", "overSampling", "predictDecisionTree", "options.filter", "majorityClass", "randomCombination", "randomWhichMax", "vector2matrix", "which.is.na", "which.is.factor", "factor2vector", "outputPerturbationSampling", "rmNA", "count.factor", "find.idx", "genericOutput", "fillVariablesNames", "is.wholenumber", "rm.tempdir", "setManyDatasets", "onlineCombineRUF", "mergeLists", "classifyMatrixCPP", "L2DistCPP", "checkUniqueObsCPP", "crossEntropyCPP", "giniCPP", "L2InformationGainCPP", "entropyInformationGainCPP", "runifMatrixCPP", "NAfactor2matrix", "factor2matrix", "as.true.matrix", "NAFeatures", "NATreatment", "rmInf", "rm.InAList")
		if (nrow(X) < 2001)
		{
			newX <- foreach( i= 1:nFeatures, .options.smp = smpopts, .inorder = FALSE, .combine = cbind, 
			.multicombine = TRUE, .export = export) %dopar%
			{
				if (X.factors[processedFeatures[i]] == 1)
				{
					rufObject <- randomUniformForest.default(X[-idx[[i]],-processedFeatures[i]], 
								Y = as.factor(X[-idx[[i]], processedFeatures[i]]), OOB = FALSE, importance = FALSE, ntree = ntree, mtry = mtry, nodesize = nodesize, threads = 1, categoricalvariablesidx = categoricalvariablesidx)
							
					X[idx[[i]], processedFeatures[i]] <- as.numeric(as.vector(predict.randomUniformForest(rufObject, 
						X[idx[[i]], -processedFeatures[i]])))
				}
				else
				{
					rufObject <- randomUniformForest.default(X[-idx[[i]],-processedFeatures[i]], Y = X[-idx[[i]], processedFeatures[i]], 
								OOB = FALSE, importance = FALSE, ntree = ntree, mtry = mtry, nodesize = nodesize, threads = 1,
								categoricalvariablesidx =  categoricalvariablesidx)
						
					X[idx[[i]], processedFeatures[i]] <- predict.randomUniformForest(rufObject, 
					X[idx[[i]], -processedFeatures[i]])	
				}							
				if (mean(is.wholenumber(rmNA(X[-idx[[i]], processedFeatures[i]]))) == 1) {  round(X[,processedFeatures[i]]) }
				else {  X[,processedFeatures[i]]  }	
			}			
			X[,processedFeatures] = newX
			stopCluster(Cl)				
		}
		else
		{
			newX <- foreach(i = 1:nFeatures, .options.smp = smpopts, .inorder = FALSE, .combine = cbind, 
			.multicombine = TRUE, .export = export) %dopar%
			{
				if (X.factors[processedFeatures[i]] == 1)
				{
					rufObject <- rUniformForest.big(X[-idx[[i]],-processedFeatures[i]], Y = as.factor(X[-idx[[i]], 
						processedFeatures[i]]), nforest = max(1, floor(nrow(X[-idx[[i]],])/2000)), replacement = TRUE, randomCut = TRUE, OOB = FALSE, importance = FALSE, ntree = ntree, mtry = mtry, 
						nodesize = nodesize, threads = 1, categoricalvariablesidx =  categoricalvariablesidx)
							
					X[idx[[i]], processedFeatures[i]] <- as.numeric(as.vector(predict.randomUniformForest(rufObject, 
						X[idx[[i]], -processedFeatures[i]])))
				}
				else
				{
					rufObject <- rUniformForest.big(X[-idx[[i]],-processedFeatures[i]], Y = X[-idx[[i]], processedFeatures[i]], 
						nforest = max(1, floor(nrow(X[-idx[[i]],])/2000)), replacement = TRUE, randomCut = TRUE, OOB = FALSE, 
						importance = FALSE, ntree = ntree, mtry = mtry, nodesize = nodesize, threads = 1, 
						categoricalvariablesidx = categoricalvariablesidx)
						
					X[idx[[i]], processedFeatures[i]] <- predict.randomUniformForest(rufObject, X[idx[[i]], -processedFeatures[i]])	
				}
							
				if (mean(is.wholenumber(rmNA(X[-idx[[i]], processedFeatures[i]]))) == 1) 
				{  round(X[,processedFeatures[i]]) }
				else  {  X[,processedFeatures[i]]  }	
			}
			X[,processedFeatures] = newX
			stopCluster(Cl)
		}			
		if (sum(X.factors) != 0)
		{
			factorsIdx = which(X.factors != 0)
			X = as.true.matrix(X)
			X = as.data.frame(X)
			for (j in 1:length(factorsIdx))
			{	
				k = factorsIdx[j]
				X[,k] = as.factor(X[,k])
				Xlevels = as.numeric(names(table(X[,k])))
				asFactorTrueX_k = as.factor(trueX[,k])
				oldLevels = levels(asFactorTrueX_k)
				for (i in 1:length(NAgrep))
				{
					levelToRemove = which(oldLevels == NAgrep[i])
					if (length(levelToRemove) > 0) 
					{ 
						if (NAgrep[i] == "")
						{	levels(X[,k]) = c("virtualClass", oldLevels[-levelToRemove])	}
						else
						{  levels(X[,k]) = oldLevels[-levelToRemove] }
					}
					else 
					{ 
						checkLevel = sample(Xlevels,1)
						if ((!is.numeric(trueX[,k])) & (checkLevel > 0) & (is.wholenumber(checkLevel)))
						{ levels(X[,k]) = oldLevels[Xlevels] }
					}
				}
			}				
		}
		for (j in 1:ncol(X))
		{
			classTrueX = class(trueX[,j])
			if ((classTrueX) == "character") X[,j] = as.character(X[,j])
		}	
		return(X)
	}
}

rufImpute = fillNA2.randomUniformForest

# unsupervised Learning = unsupervised2supervised + rUF model + proximity + MDS/Spectral (+ gap) + kMeans/hClust 
unsupervised2supervised <- function(X, 
method = c("uniform univariate sampling", "uniform multivariate sampling"),
seed = 2014,
conditionalTo = NULL,
samplingFromGaussian = FALSE,  
bootstrap = FALSE)
{
	np = dim(X); n = np[1]; p = np[2]
	flag = FALSE
	syntheticTrainLabels1 = rep(0,n)
	if (!is.null(rownames(X))) { rowNamesX = rownames(X); flag = TRUE }
	if (is.data.frame(X))
	{	
		X = NAfactor2matrix(X, toGrep = "anythingParticular") 
		cat("X is a data frame and has been converted to a matrix.\n")
	}
	
	XX = matrix(data = NA, nrow = n, ncol = p)
	
	set.seed(seed)
	if (method[1] == "uniform univariate sampling")
	{	
		if (is.null(conditionalTo))
		{	XX <- apply(X, 2, function(Z) sample(Z, n, replace = bootstrap))  }
		else
		{
			if (length(conditionalTo) == nrow(X)) {	 ZZ = round(as.numeric(conditionalTo),0) }
			else {	ZZ = round(as.numeric(X[,conditionalTo[1]]),0) }
			ZZtable = as.numeric(names(table(ZZ)))
			XX = X
			for (i in 1:length(ZZtable))
			{
				idx = which(ZZ == ZZtable[i])
				if (samplingFromGaussian)
				{  XX[idx,] <- apply(X[idx,], 2, function(Z) rnorm(length(idx), mean(Z), sd(Z)))  }
				else				
				{	XX[idx,] <- apply(X[idx,], 2, function(Z) sample(Z, length(idx), replace = bootstrap))	}
			}
		}
		
	}
	else
	{
		if (is.null(conditionalTo))
		{
			for (j in 1:p)
			{
				XX[,j] = sample(X, n, replace = bootstrap)
				
				outOfRange = which( (XX[,j] > max(X[,j])) | (XX[,j] < min(X[,j])) )
				while (length(outOfRange) > 0)
				{	
					XX[outOfRange, j] = sample(X, length(outOfRange), replace = bootstrap)	
					outOfRange = which( (XX[,j] > max(X[,j])) | (XX[,j] < min(X[,j])) )
				}
			}
		}
		else
		{
			if (length(conditionalTo) == nrow(X)) {	 ZZ = round(as.numeric(conditionalTo),0) }
			else {	ZZ = round(as.numeric(X[,conditionalTo[1]]),0) }
			ZZtable = as.numeric(names(table(ZZ)))
			XX = X
			meanX = mean(X)
			sdX = sd(X)
			for (i in 1:length(ZZtable))
			{
				idx = which(ZZ == ZZtable[i])
				for (j in 1:p)
				{
					if (samplingFromGaussian)
					{  XX[idx,j] = rnorm(length(idx), meanX, sdX)	}
					else 
					{ 	XX[idx,j] = sample(X, length(idx), replace = bootstrap)	}
					
					outOfRange = which( (XX[idx,j] > max(X[,j])) | (XX[idx,j] < min(X[,j])) )
					while (length(outOfRange) > 0)
					{	
						XX[idx[outOfRange], j] = sample(X, length(outOfRange), replace = bootstrap)	
						outOfRange = which( (XX[idx[outOfRange],j] > max(X[,j])) | 
						(XX[idx[outOfRange],j] < min(X[,j])) )
					}
				}
			}
		}
	}
	newX = rbind(X, XX)
	if (flag) { rownames(newX) = c(rowNamesX, rowNamesX) }
	syntheticTrainLabels2 = rep(1, n)	
	Y = c(syntheticTrainLabels1, syntheticTrainLabels2)
	
	return(list(X = newX, Y = Y))	
}

proximitiesMatrix <- function(object, fullMatrix = TRUE, Xtest = NULL, predObject = NULL, sparseProximities = TRUE, pthreads = "auto")
{ 
	# proximity matrix of cases : the bigger the value, the closest is proximity of i an j cases
	object = filter.object(object)

	if (is.null(object$predictionObject))
	{
		if (is.null(predObject)) 
		{ 
			if (is.null(Xtest))	{ stop("please provide test data.\n") }
			else
			{
				predObject <- predict.randomUniformForest(object, Xtest, type = "all")
				votesData = predObject$votes.data
			}
		}
		else 
		{  
			if (is.null(predObject$majority.vote)) 
			{ stop("Please provide full prediction object (type = 'all') when calling predict() function") }
			else 
			{	votesData = predObject$votes.data 	}
		}				
	}		
	else
	{	votesData = object$predictionObject$votes.data	}
		
	n = nrow(votesData[[1]])  
	B = length(votesData)
	
	if ((n < 500) & is.character(pthreads)) { pthreads = 1 }
	else	
	{
		#require(parallel)	
		max_threads = detectCores()					
		if (pthreads == "auto")
		{	
			if (max_threads == 2) { pthreads = max_threads }
			else {	pthreads  = max(1, max_threads - 1)  }
		}
		else
		{
			if (max_threads < pthreads) 
			{  cat("Warning : number of threads is higher than logical threads in this computer.\n") }
		}			
		
		#require(doParallel)	
		Cl = makePSOCKcluster(pthreads, type = "SOCK")
		registerDoParallel(Cl)
	}
			
	if (fullMatrix)
	{	
		operatorLimit = 3*800*n^2/(10000^2)
		if ( memory.limit() < operatorLimit)
		{ 
			cat("Proximity matrix is likely to be too big. Model will compute a 'n x B' matrix with B = number of trees\n")
			proxMatrix = matrix(0L, B, n)
			fullMatrix = FALSE
		}
		else
		{	proxMatrix = matrix(0L, n, n) }
		
		proxMat = proxMatrix
		if (!sparseProximities)
		{
			if (pthreads == 1)
			{
				for (i in 1:n)
				{		
					for (b in 1:B)
					{
						common.idx = which(votesData[[b]][i,3] == votesData[[b]][,3])
						proxMatrix[i, common.idx] = 1L + proxMatrix[i,common.idx] 
					}
				}			
			}
			else
			{
				proxMatrix <- foreach(i = 1:n, .combine = rbind, .multicombine = FALSE) %dopar%
				{		
					for (b in 1:B)
					{
						common.idx = which(votesData[[b]][i,3] == votesData[[b]][,3])
						proxMat[i, common.idx] = 1L + proxMat[i, common.idx] 
					}
					proxMat[i,] 
				}
				#proxMatrix =proxMat +
				stopCluster(Cl)
			}
		}
		else
		{
			if (pthreads == 1)
			{
				for (i in 1:n)
				{		
					for (b in 1:B)
					{
						common.idx = which(votesData[[b]][i,3] == votesData[[b]][,3])
						proxMatrix[i, common.idx] = 1L + proxMat[i,common.idx] 
					}
				}			
			}
			else
			{
				proxMatrix <- foreach(i = 1:n, .combine = rbind, .multicombine = TRUE) %dopar%
				{		
					for (b in 1:B)
					{
						common.idx = which(votesData[[b]][i,3] == votesData[[b]][,3])
						proxMat[i, common.idx] = 1L + proxMatrix[i, common.idx] 
					}
					proxMat[i,] 
				}
				#proxMatrix =proxMat +
				stopCluster(Cl)
			}
		}
	}
	else
	{	
		proxMatrix = matrix(0L, n, B)
		proxMat = proxMatrix
		
		if (pthreads == 1)
		{
			for (i in 1:n)
			{		
				for (b in 1:B)
				{
					common.idx = which((votesData[[b]][i,3] == votesData[[b]][,3]))
					proxMat[common.idx, b] = 1L + proxMatrix[common.idx,b] 
				}
				proxMat[i,]
			}
		}
		else
		{
			proxMatrix <- foreach(i = 1:n, .combine = rbind, .multicombine = TRUE) %dopar%
			{		
				for (b in 1:B)
				{
					common.idx = which((votesData[[b]][i,3] == votesData[[b]][,3]))
					proxMat[common.idx, b] = 1L + proxMatrix[common.idx,b] 
				}
				proxMat[i,]
			}
			stopCluster(Cl)
		}		
	}
	varNames = NULL
	if (fullMatrix) { nn = n }
	else { nn = B }
	for (i in 1:nn) { varNames = c(varNames,paste("C", i, sep="")) }
              	
	colnames(proxMatrix) = varNames
	rownames(proxMatrix) = if (!is.null(rownames(Xtest))) { rownames(Xtest) } else { 1:n }
	
	if (sparseProximities) { return(proxMatrix) } else { return(proxMatrix/B) }
	#return(proxMatrix/B)
}
	
kBiggestProximities <- function(proximitiesMatrix, k) t(apply(proximitiesMatrix, 1, function(Z) sort(Z, decreasing = TRUE)[1:(k+1)][-1] ))

MDSscale <- function(proximitiesMatrix, metric = c("metricMDS", "nonMetricMDS"), dimension = 2,  distance = TRUE, plotting = TRUE, seed = 2014)
{
	set.seed(seed)
	if (metric[1] == "metricMDS")  
	{ 
		if (distance)  { d = stats::dist(1 - proximitiesMatrix) }
		else 
		{ 
			if (nrow(proximitiesMatrix) == ncol(proximitiesMatrix)) {	d = 1 - proximitiesMatrix  }
			else { d = stats::dist(1 - proximitiesMatrix); cat("Matrix is not a square matrix. Distance is used.\n") }
		}
		fit = stats::cmdscale(d, eig = TRUE, k = dimension) 
	}	
	else
	{	
		eps = 1e-6
		if (distance)  { d = eps + stats::dist(1 - proximitiesMatrix) }
		else 
		{ 
			if (nrow(proximitiesMatrix) == ncol(proximitiesMatrix)) {	d = 1 - (proximitiesMatrix - eps)  }
			else { d = eps + stats::dist(1 - proximitiesMatrix); cat("Matrix is not a square matrix. Distance is used.\n") }
		}
		
		#require(MASS)
		fit = MASS::isoMDS(d, k = dimension)
	}
		
	if (plotting)
	{ 
		x <- fit$points[,1]
		y <- fit$points[,2]
		plot(x, y, xlab = "Coordinate 1", ylab = "Coordinate 2", main = metric[1], type = "p", lwd = 1, pch = 20)
	}
	
	return(fit)
}

specClust <- function(proxMat, k = floor(sqrt(ncol(proxMat))))
{
	A = proxMat
	diag(A) = 0
	np = dim(A)
	n = np[1]
	p = np[2]
	D = matrix(0, n, p)
	diag(D) = (rowSums(A))^(-0.5)
	L = (D%*%A%*%D)
	
	X = eigen(L)$vectors[,1:k]
	Xcol = (colSums(X^2))^0.5

	return(X/Xcol)
}

gap.stats <- function(fit, B = 100, maxClusters = 5, seed = 2014)
{
	set.seed(seed)
	#require(cluster)
	return(cluster::clusGap(fit, FUN = kmeans, K.max = max(2, maxClusters), B = B))	
}

kMeans <- function(fit, gapStatFit, maxIters = 10, plotting = TRUE, algorithm = NULL, k = NULL, 
reduceClusters = FALSE, seed = 2014)
{
	set.seed(seed)
	if (is.null(maxIters)) { maxIters = 10 }
	if (is.null(algorithm)) { algorithm = "Hartigan-Wong" }
	
	if (is.null(k))
	{
		k = cluster::maxSE(gapStatFit$Tab[, "gap"], gapStatFit$Tab[, "SE.sim"])
		if (k == 1) { cat("Only one cluster has been found.\nNumber of clusters has been set to 2.\n") ; k = 2 }
		if (k == length(gapStatFit$Tab[, "gap"])) 
		{  cat("\nNumber of clusters found is equal to the maximum number of clusters allowed.\n") }
	}
	fit.kMeans <- stats::kmeans(fit, k, iter.max = maxIters, algorithm = algorithm)
	
	if ( (k > 2) & reduceClusters)
	{
		clusterSizes = table(fit.kMeans$cluster)
		smallSizes = which(clusterSizes < 0.05*length(fit[,1]))
		nSmallSizes = length(smallSizes)
		if (nSmallSizes > 0)
		{ 
			smallSizes = sort(smallSizes, decreasing = TRUE)
			for (i in 1:nSmallSizes)
			{	
				fit.kMeans$cluster[which(fit.kMeans$cluster == smallSizes[i])] =  as.numeric(names(clusterSizes))[max(1,smallSizes-1)]
			}
		}		
	}
	
	if (plotting)
	{
		x = fit[,1]
		y = fit[,2]
		plot(x, y, xlab = "Coordinate 1", ylab = "Coordinate 2", main = "Multi-dimensional scaling", type = "p", lwd = 1, pch = 20)
		for (i in 1:k)
		{  
			idx = which(fit.kMeans$cluster == i)
			points(x[idx], y[idx], type = "p", lwd = 1, pch = 20, col = i)
			#legend(max(x) - 1.25, max(y) + 0.35 + 0.5*(1-i), paste("Cluster ", i), fill= i, border =NA, bty="n", col = i)
		}
	}
	
	return(fit.kMeans)
}

hClust <- function(proximities, method = NULL, plotting = TRUE, k = NULL, reduceClusters = FALSE, seed = 2014)
{
	set.seed(seed)
	if (is.null(method)) { method = "complete" }
	d <- stats::dist(1 - proximities)
	XClust <- stats::hclust(d, method = method,  members = NULL)
	
	if (is.null(k))
	{
		diffHeightsIdx = which.max(diff(XClust$height))
		height = (XClust$height[diffHeightsIdx] +  XClust$height[diffHeightsIdx-1])/2
		values <- stats::cutree(XClust, h = height)
		nbClusters = length(unique(values))
	}
	else
	{
		nbClusters = k
		values <- stats::cutree(XClust, k = nbClusters)
		height = NULL
	}
	
	if ( (nbClusters > 2) & reduceClusters)
	{
		clusterSizes = table(values)
		smallSizes = which(clusterSizes < 0.05*length(values))
		nSmallSizes = length(smallSizes)
		if (nSmallSizes > 0)
		{ 
			smallSizes = sort(smallSizes, decreasing = TRUE)
			for (i in 1:nSmallSizes)
			{	
				values[which(values == smallSizes[i])] = as.numeric(names(clusterSizes))[max(1,smallSizes - 1)]
			}
		}		
	}
	
	if (plotting)
	{	
		showLabels = if (nrow(proximities) < 300) { TRUE } else { FALSE }
		plot(XClust, xlab = "Observations", labels = showLabels)
				# A2Rplot(XClust, k = nbClusters, boxes = TRUE, show.labels = FALSE, main = "Dendrogram")
		if (!is.null(height)) { abline(h = height , col='red') }
	}
	
    return(list(object = XClust, cluster = values))
}

observationsImportance <- function(X, importanceObject)
{
	X = factor2matrix(X)
	object = importanceObject$localVariableImportance$obsVariableImportance[,-1]
	idx = grep("localVariableFrequency",colnames(object))
	object = object[,-idx]
	n = nrow(X)
	XX = matrix(NA, n, ncol(object))
	for (i in 1:nrow(X))
	{	XX[i,] = X[i, object[i,]]	}
	
	rownames(XX) = 1:n
	return(XX)
}

print.unsupervised <- function(x, ...)
{
	object = x
	rm(x)
	Z = object$unsupervisedModel$cluster
	if (!is.factor(Z))
	{
		if (!is.null(object$unsupervisedModel$clusterOutliers))
		{	
			Z = c(Z, object$unsupervisedModel$clusterOutliers)				
			Z = sortDataframe( data.frame(Z, as.numeric(names(Z))), 2)
			Z = Z[,1]
		}
	}	
	ZZ = table(Z)
	names(attributes(ZZ)$dimnames) = NULL
	
	if (object$params["endModel"] == "MDS")
	{
		x = object$MDSModel$points[,1]
		y = object$MDSModel$points[,2]
		percentExplained = (interClassesVariance(x, Z) +  interClassesVariance(y, Z))/(variance(x) + variance(y))
		
		cat("Average variance between clusters (in percent of total variance): ", round(100*percentExplained,2), "%\n", sep = "")
				
		cat("Average silhouette: ", round(mean(cluster::silhouette(as.numeric(Z), dist(object$MDSModel$points))[,3]),4), "\n",  sep = "")
		cat("Clusters size:\n")
		print(ZZ)
		cat("Clusters centres (in the MDS coordinates):\n") 
		print(round(object$unsupervisedModel$centers,4))		
	}
		
	if ((object$params["endModel"] == "MDSkMeans") | (object$params["endModel"] == "SpectralkMeans"))
	{
		cat("Average variance between clusters (in percent of total variance): ", round(100*sum(object$unsupervisedModel$betweenss)/sum(object$unsupervisedModel$totss),2), "%\n", sep = "")
				
		cat("Average silhouette: ", round(mean(cluster::silhouette(Z, dist(object$MDSModel$points))[,3]),4), "\n",
		sep = "")
		cat("Clusters size:\n")
		print(ZZ)
		cat("Clusters centres (in the MDS coordinates):\n") 
		print(round(object$unsupervisedModel$centers,4))		
	}
	
	if (object$params["endModel"] == "MDShClust")
	{
		x = object$MDSModel$points[,1]
		y = object$MDSModel$points[,2]
		
		percentExplained1 = interClassesVariance(x, Z)/variance(x)
		percentExplained2 = interClassesVariance(y, Z)/variance(y)
		percentExplained  = round(0.5*(percentExplained1 + percentExplained2), 4)				
		cat("Average variance between clusters (in percent of total variance) : ", percentExplained*100, "%\n", sep = "")
		
		cat("Average silhouette : ", round(mean(cluster::silhouette(Z, dist(object$MDSModel$points))[,3]),4), "\n",  sep = "")
		cat("Clusters size:\n")
		print(ZZ)
		cat("Clusters centers (in the MDS coordinates):\n") 
		print(round(object$unsupervisedModel$centers,4))		
	}
}

plot.unsupervised <- function(x, importanceObject = NULL, xlim = NULL, ylim = NULL, coordinates = NULL,...)
{
	object = x
	rm(x)
	if (is.null(coordinates)) 
	{ 
		coordinates = c(as.numeric( substr(object$params["coordinates"], 1, 1)), 
		as.numeric( substr(object$params["coordinates"], 2, 2)) )
		if (ncol(object$MDSModel$points) == 2) 
		{ 
			if (coordinates[1] > 1)
			{ 
				cat(paste("Coordinates", coordinates[1], "and", coordinates[2], "have been plotted, but appear as Coordinate 1 and Coordinate 2.\n", sep=" "))
			}
			coordinates = 1:2 
		}		
	}
	else
	{
		if (ncol(object$MDSModel$points) == 2) coordinates = 1:2 
	}

	x = object$MDSModel$points[,coordinates[1]]
	y = object$MDSModel$points[,coordinates[2]]
	
	offsetY = 0
	if (!is.null(importanceObject))
	{
		clusterFeatures = importanceObject$localVariableImportance$classVariableImportance
		p = ncol(clusterFeatures)
		clusterFeaturesNames = NULL
		varNames = varValues = vector(length = p)
		for (j in 1:p)
		{
			clusterFeatures = sortDataframe(clusterFeatures, j, decrease = TRUE)
			idx = which(clusterFeatures[,j] > 0.05)
			varNames[j] = rownames(clusterFeatures)[idx[1]]
			varValues[j] = clusterFeatures[idx[[1]],j]
			clusterFeaturesNames = c(clusterFeaturesNames, paste(varNames[j], ":", round(100*varValues[j],0), "%",	sep ="" ))
		}
		offsetY = -abs(diff(range(y)))*0.1
	}
		
	clusters = object$unsupervisedModel$cluster
	
	if (object$params["endModel"] != "MDS")
	{
		if (!is.null(object$unsupervisedModel$clusterOutliers))
		{	
			clusters = c(clusters, object$unsupervisedModel$clusterOutliers)				
			clusters = sortDataframe( data.frame(clusters, as.numeric(names(clusters))),2)
			clusters = clusters[,1]
		}
	}
			
	offsetX = abs(diff(range(x)))*0.3
	uniqueClusters = sort(unique(clusters))
	
	dev.new()
	if (object$params["endModel"] == "MDShClust")
	{
		dev.off()
		#percentExplained = 0.5*(percentExplained1 + percentExplained2)
		nbClusters =  length(uniqueClusters)
		
		XClust = object$unsupervisedModel$object
		diffHeightsIdx = which.max(diff(XClust$height))
		height = (XClust$height[diffHeightsIdx] +  XClust$height[diffHeightsIdx-1])/2				
		showLabels = if (nrow(object$proximityMatrix) < 300) { TRUE } else { FALSE }
		plot(XClust, labels = showLabels)
		abline(h = height , col='red') 
		
		# A2Rplot(object$unsupervisedModel$object, k = nbClusters, boxes = TRUE, show.labels = FALSE, 
		# main = "Dendrogram", col.down = nbClusters:1)
		
		cat("Two graphics have been plotted. Please slide the window to see the second one.\n")
		dev.new()		
	}
	
	if (is.null(xlim)) { xlim = c(min(x) - offsetX, max(x) + offsetX) }
	if (is.null(ylim)) { ylim = c(min(y) + offsetY, max(y)) }
	
	if (length(grep("MDS", object$params["endModel"])) == 1)  
	{ plotTitle =  "Multidimensional scaling and Clusters representation" }
	else 
	{  plotTitle =  "Spectral decomposition and Clusters representation" }
	
	plot(x, y, xlab = paste("Coordinate", coordinates[1]), ylab = paste("Coordinate", coordinates[2]), 
	main = plotTitle, type = "p", lwd = 1, pch = 20, xlim = xlim, ylim = ylim)
	nbCusters = length(uniqueClusters)
	for (i in 1:nbCusters)
	{  
		idx = which(clusters == uniqueClusters[i])
		points(x[idx], y[idx], type = "p", lwd = 1, pch = 20, col = i)
		#if (object$params["endModel"] != "MDShClust")
		{	
			points(object$unsupervisedModel$centers[i,coordinates[1]], 
			object$unsupervisedModel$centers[i,coordinates[2]], lwd = 1, cex = 1.5, pch = 8, col = i)	
		}
	}
		
	if (!is.null(importanceObject))
	{ 	
		clusterNames = rm.string(colnames(clusterFeatures), "Class ")
		legend("topright", inset = .01, clusterNames, fill = 1:i, horiz = FALSE, border = NA, bty ="n")
		legend("bottomleft", cex = 0.7, as.character(clusterFeaturesNames), fill = 1:i, horiz = FALSE, border = NA, bty ="n")		
		print(clusterFeatures)		
	}
	else
	{ legend("topright", inset = .01, as.character(uniqueClusters), fill = 1:i, horiz = FALSE, border = NA, bty ="n")}
	
}

modifyClusters <- function(object, decreaseBy = NULL, increaseBy = NULL, seed = 2014)
{
	params = object$modelParams
	endModel = object$params["endModel"]
	nbClusters = object$nbClusters
	
	if (!is.null(decreaseBy))
	{	k = nbClusters - decreaseBy }
	
	if (!is.null(increaseBy))
	{	k = nbClusters + increaseBy }
	
	if (is.null(decreaseBy) & (is.null(increaseBy)))
	{	stop("Please increase or decrease number of clusters")  }
	
	if (k <= 1) 
	{ 
		cat("Only one cluster will remain. Model has not been modified.\n") 
		return(object)
	}
	else
	{
		endModelMetric = if (object$params["endModelMetric"] == "NULL") { NULL } 
		else { object$params["endModelMetric"] }
		Z = object$MDSModel$points
		
		if ( (endModel == "MDSkMeans") | (endModel == "SpectralkMeans") )
		{
			maxIters = if (object$params["maxIters"] == "NULL") { 10 } else { as.numeric(object$params["maxIters"]) }
			X.model <- kMeans(Z, NULL, k = k, maxIters = maxIters, algorithm = endModelMetric, plotting = FALSE, reduceClusters = FALSE, seed = seed)
		}
		
		if (endModel == "MDShClust")
		{	
			X.model <- hClust(Z, method = endModelMetric, plotting = FALSE, k = k, reduceClusters = FALSE, seed = seed)
			centers = matrix(NA, k, ncol(Z))
			for (i in 1:k)
			{
				idx = which(X.model$cluster == i)
				centers[i,] = colMeans(Z[idx,, drop = FALSE])
			}
			X.model$centers = centers
		}
		object$unsupervisedModel = NULL
		object$unsupervisedModel = X.model  
		object$nbClusters = k
		
		return(object)
	}
}

splitClusters <- function(object, whichOnes, seed = 2014, ...) 
{
	cat("Note that outliers are not currently taken into account.\n")
	whichOnes = sort(whichOnes)
	
	endModel = object$params["endModel"]
	if (exists("endModelMetric")) { endModelMetric = endModelMetric }
	else
	{
		if (object$params["endModelMetric"] == "NULL") { endModelMetric = NULL } 
		else { endModelMetric = object$params["endModelMetric"] }	
	}
	Z = object$MDSModel$points
	nClusters = length(unique(object$unsupervisedModel$cluster))
		
	for (i in 1:length(whichOnes))
	{
		idx = which(object$unsupervisedModel$cluster == whichOnes[i])	
		if (length(idx) == 0) { stop(paste("No elements available for cluster", i,".\n")) }
		
		if ( (endModel == "MDSkMeans") | (endModel == "SpectralkMeans") )
		{
			if (exists("maxIters")) { maxIters = maxIters }
			else
			{
				if (object$params["maxIters"] == "NULL") { maxIters = 10 } 
				else { maxIters = as.numeric(object$params["maxIters"]) }
			}
			X.model <- kMeans(Z[idx,], NULL, k = 2, maxIters = maxIters, algorithm = endModelMetric, plotting = FALSE, reduceClusters = FALSE, seed = seed)	
			object$unsupervisedModel$betweenss = 0.5*(object$unsupervisedModel$betweenss + X.model$betweenss)
			object$unsupervisedModel$totss = 0.5*(object$unsupervisedModel$totss + X.model$totss)
		}
		
		if (endModel == "MDShClust")
		{	
			X.model <- hClust(Z, method = endModelMetric, plotting = FALSE, k = 2, reduceClusters = FALSE, seed = seed)
			centers = matrix(NA, 2, ncol(Z))
			for (j in 1:2)
			{
				idx2 = which(X.model$cluster == j)
				centers[j,] = colMeans(Z[idx2, , drop = FALSE])
			}
			X.model$centers = centers
		}
		  
		idx3 = which(X.model$cluster == 2)
		object$unsupervisedModel$cluster[idx][idx3] = nClusters + i
	}
	
	nClusters = length(unique(object$unsupervisedModel$cluster))
	centers = matrix(NA, nClusters, ncol(Z))
	for (cl in 1:nClusters)
	{
		idx = which(object$unsupervisedModel$cluster == cl)
		centers[cl,] = colMeans(Z[idx, , drop = FALSE])
	}
	
	object$unsupervisedModel$centers = centers
	rownames(object$unsupervisedModel$centers) = 1:nClusters
		
	return(object)
}

mergeClusters <- function(object, whichOnes)
{
	whichOnes = sort(whichOnes)
	
	idx1 = which(object$unsupervisedModel$cluster == whichOnes[1])	
	object$unsupervisedModel$cluster[idx1] = whichOnes[2]
	
	flag = FALSE
	if (!is.null(object$unsupervisedModel$clusterOutliers))
	{  
		idx1 = which(object$unsupervisedModel$clusterOutliers == whichOnes[1])
		if (length(idx1) > 0) { object$unsupervisedModel$clusterOutliers[idx1] = whichOnes[2] }
		flag = TRUE
	}
		
	idx2 = which(object$unsupervisedModel$cluster == whichOnes[2])
	object$unsupervisedModel$centers[whichOnes[2],] = t(colMeans(object$MDSModel$points[idx2,]))
	object$unsupervisedModel$centers = object$unsupervisedModel$centers[-whichOnes[1],]
	
	idx = sort(unique(object$unsupervisedModel$cluster))
	for (i in 1:length(idx))
	{
		newIdx = which(object$unsupervisedModel$cluster == idx[i])
		object$unsupervisedModel$cluster[newIdx] = i
		if (flag)
		{  
			newOutliersIdx = which(object$unsupervisedModel$clusterOutliers == idx[i])
			if (length(newOutliersIdx) > 0)  { object$unsupervisedModel$clusterOutliers[newOutliersIdx] = i }
		}
	}
	rownames(object$unsupervisedModel$centers) = 1:i
	cat(paste("\ncluster ", idx, " has been changed to ", 1:i, sep=""))
	cat("\n")
	
	return(object)
}

clusteringObservations <- function(object, X, OOB = TRUE, predObject = NULL, importanceObject = NULL, 
baseModel = c("proximity", "proximityThenDistance", "importanceThenDistance"),
MDSmetric = c("metricMDS", "nonMetricMDS"),
endModel = c("MDS", "MDSkMeans", "MDShClust", "SpectralkMeans"),
...)
{
	clusterObject <- unsupervised.randomUniformForest(object, Xtest = X, importanceObject = importanceObject, predObject = predObject, baseModel = baseModel[1], MDSmetric = MDSmetric[1], endModel = endModel[1], 
	outliersFilter = FALSE, OOB = OOB, ...)
	plot(clusterObject, importanceObject = importanceObject, ...)
	return(clusterObject)
}

as.supervised <- function(object, X, ...)
{
	if (class(object) != "unsupervised")
	{	stop("Please provide an unsupervised randomUniformForest object.\n") }
	
	n = nrow(X)
	clusters = object$unsupervisedModel$cluster
	if (!is.null(object$unsupervisedModel$clusterOutliers))
	{	clusters = c(clusters, object$unsupervisedModel$clusterOutliers)	}
	if (is.null(names(clusters))) { names(clusters) = 1:n }
	clusters  = sortDataframe( data.frame(clusters, as.numeric(names(clusters ))),2)
	Y = as.factor(clusters[,1])
	
	if (n < 10000) { rUFObject <- randomUniformForest(X, Y, ...)	}
	else { rUFObject <- rUniformForest.big(X, Y, nforest = min(2, floor(n/2000)), randomCut = TRUE, ...)  }
	
	return(rUFObject)	
}

update.unsupervised <- function(object, X = NULL, oldData = NULL, mapAndReduce = FALSE, updateModel = FALSE, ...)
{
	if(is.null(X)) { stop("'X' is null. Please provide new data.") }
	
	endModel = object$params["endModel"]
	clusters = as.numeric(object$params["clusters"])
	maxIters = if (object$params["maxIters"] == "NULL") NULL else as.numeric(object$params["maxIters"])
	endModelMetric = if (object$params["endModelMetric"] == "NULL") NULL else object$params["endModelMetric"]
	reduceClusters = if (object$params["reduceClusters"] == "FALSE") FALSE else TRUE
	seed = as.numeric(object$params["seed"])
	
	if ( (endModel != "MDShClust") & (endModel != "MDSkMeans") & (endModel != "SpectralkMeans") )
	{ 
		stop("Update is only available for models which were using endModel = 'MDShClust', endModel = 'MDSkMeans' or endModel = 'SpectralkMeans' options") 
	}
	
	p = ncol(object$MDS$points)
	updateLearningMDS = predictedMDS = learningMDS = vector('list', p)
	if (is.null(object$largeDataLearningModel))
	{
		if (is.null(oldData))
		{ stop("Former data are needed to (learn MDS points and) update model.\n") }
		else
		{
			if (mapAndReduce)
			{	
				for (j in 1:p)
				{ learningMDS[[j]] = rUniformForest.big(oldData, object$MDSModel$points[,j], randomCut = TRUE, ...) } 
			}
			else
			{
				for (j in 1:p)
				{ learningMDS[[j]] = randomUniformForest(oldData, object$MDSModel$points[,j], ...)  }
			}
		}
	}
	else
	{
		for (j in 1:p)
		{ learningMDS[[j]] = object$largeDataLearningModel[[j]] }
	}
	
	formerMDSPoints = object$MDSModel$points
	for (j in 1:p)
	{  predictedMDS[[j]] = predict(learningMDS[[j]], X) }
	
	newMDSPoints = do.call(cbind, predictedMDS)		
	Z = rbind(formerMDSPoints, newMDSPoints)
	
	if (updateModel)
	{
		if (mapAndReduce)
		{	
			for (j in 1:p)
			{  updateLearningMDS[[j]] = rUniformForest.big(X, predictedMDS[[j]], randomCut = TRUE, ...) }
		}
		else
		{
			for (j in 1:p)	{ updateLearningMDS[[j]] = randomUniformForest(X, predictedMDS[[j]], ...) } 
		}
		
		for (j in 1:p) 	{ learningMDS[[j]] = rUniformForest.combine(learningMDS[[j]], updateLearningMDS[[j]]) }
	}
	
	if ( (endModel[1] == "MDSkMeans") | (endModel[1] == "SpectralkMeans") )
	{
		X.model <- kMeans(Z, NULL, k = clusters, maxIters = maxIters, plotting = FALSE, algorithm = endModelMetric, reduceClusters = reduceClusters, seed = seed)			
	}
			
	if (endModel[1] == "MDShClust")  
	{	
		X.model <- hClust(Z, method = endModelMetric, plotting = FALSE, k = clusters, 
		reduceClusters = reduceClusters,  seed = seed)
		centers = matrix(NA, clusters, ncol(Z))
		for (i in 1:clusters)
		{
			idx = which(X.model$cluster == i)
			centers[i,] = colMeans(Z[idx,, drop = FALSE])
		}
		X.model$centers = centers
	}
	
	object$X.scale$points = Z
	largeDataLearningModel = learningMDS
	unsupervisedObject = list(proximityMatrix = object$proxMat, MDSModel = object$X.scale, 
		unsupervisedModel = X.model, largeDataLearningModel = largeDataLearningModel, gapStatistics = NULL, 
		rUFObject = object$rUF.model, nbClusters = clusters, params = object$params) 	
	
	class(unsupervisedObject) <- "unsupervised"
	
	unsupervisedObject
}

combineUnsupervised <- function(...)
{
	object <- list(...)
	n = length(object)
	
	i = 1
	endModel = object[[i]]$params["endModel"]
	clusters = as.numeric(object[[i]]$params["clusters"])
	maxIters = if (object[[i]]$params["maxIters"] == "NULL") NULL else as.numeric(object[[i]]$params["maxIters"])
	endModelMetric = if (object[[i]]$params["endModelMetric"] == "NULL") NULL else object[[i]]$params["endModelMetric"]
	reduceClusters = if (object[[i]]$params["reduceClusters"] == "FALSE") FALSE else TRUE
	seed = as.numeric(object[[i]]$params["seed"])
	
	if ( (endModel != "MDShClust") & (endModel != "MDSkMeans") & (endModel != "SpectralkMeans") )
	{ 
		stop("Combine is only available for models which were using endModel = 'MDShClust', endModel = 'MDSkMeans' or endModel = 'SpectralkMeans' options") 
	}
		
	Z = NULL
	for (i in 1:n)
	{   Z = rbind(Z, object[[i]]$MDSModel$points)  }
		
	if ( (endModel[1] == "MDSkMeans") | (endModel != "SpectralkMeans"))
	{
		X.model <- kMeans(Z, NULL, k = clusters, maxIters = maxIters, plotting = FALSE, algorithm = endModelMetric, reduceClusters = reduceClusters, seed = seed)			
	}
			
	if (endModel[1] == "MDShClust")  
	{	
		X.model <- hClust(Z, method = endModelMetric, plotting = FALSE, k = clusters, 
		reduceClusters = reduceClusters,  seed = seed)
		X.gapstat = NULL	
	}
	
	X.scale = list()
	X.scale$points = Z
		
	unsupervisedObject = list(proximityMatrix = NULL, MDSModel = X.scale, 
		unsupervisedModel = X.model, largeDataLearningModel = NULL, gapStatistics = NULL, 
		rUFObject = NULL, nbClusters = clusters, params = object[[1]]$params) 	
	
	class(unsupervisedObject) <- "unsupervised"
	
	unsupervisedObject
}

scalingMDS <- function(object)
{
	if (is.null(object$MDSModel$standardizedPoints))
	{	object$MDSModel$points = standardize(object$MDSModel$points) }
	else
	{   object$MDSModel$points = object$MDSModel$standardizedPoints 	}
	
	return(object)
}	

updateCombined.unsupervised <- function(object, X, mapAndReduce = FALSE, ...)
{
	if (!is.null(object$largeDataLearningModel))
	{	stop("The function can only rebuild a learning model of MDS (or spectral) points for formerly combined objects.\n") }
	
	if (nrow(X) !=  nrow(object$MDSModel$points))
	{ 	stop("MDS (or spectral) points and data do not have the same number of rows.\n") }		
	
	p = ncol(object$MDSModel$points)
	learningMDS = vector('list', p)
	if (mapAndReduce)
	{	
		for (j in 1:p)
		{	learningMDS[[j]] = rUniformForest.big(X, object$MDSModel$points[,j], randomCut = TRUE, ...) }
	}
	else
	{
		for (j in 1:p)
		{	learningMDS[[j]] = randomUniformForest(X, object$MDSModel$points[,j], ...) }
	}
	
	unsupervisedObject = list(proximityMatrix = NULL, MDSModel = object$MDSModel, 
		unsupervisedModel = object$unsupervisedModel, largeDataLearningModel = learningMDS, 
		gapStatistics = object$gapStatistics, rUFObject = object$rUFObject, nbClusters = object$nbClusters, 
		params = object$params) 	
	
	class(unsupervisedObject) <- "unsupervised"
	
	unsupervisedObject
}

clusterAnalysis <- function(object, X, components = 2, maxFeatures = 2, clusteredObject = NULL, categorical = NULL, 
OOB = FALSE)
{
	flagList = FALSE
	if (!is.null(clusteredObject)) 
	{ 
		if (is.list(clusteredObject)) { flagList = TRUE }
	}
	else { flagNULL = FALSE   }
		
	
	p = ncol(X)
	p.new = floor((ncol(object$localVariableImportance$obsVariableImportance)-1)/2)
	if (p.new < components) 
	{ 
		cat("Maximal number of components is set to 2. Increase 'maxInteractions' option\nin 'importance()' function in order to assign more components.\n")
		components = p.new 
	}
	n = nrow(X)
	varIdx = c(1:(components+1), (p.new+2):(p.new + components + 1))
	   
	if (class(clusteredObject)[1]  == "unsupervised") { clusterName = "Clusters" } 
	else { clusterName = "Class" }
	
	featuresAndObs = as.data.frame(object$localVariableImportance$obsVariableImportance)
	frequencyFeaturesIdx = grep("Frequency", colnames(featuresAndObs))
	featuresNames = apply(featuresAndObs[,-c(1,frequencyFeaturesIdx)], 2, function(Z) colnames(X)[Z])
	featuresAndObs[,-c(1,frequencyFeaturesIdx)] = featuresNames

	groupedAnalysis = aggregate(featuresAndObs[,2:min((components+1), p.new+1)], list(featuresAndObs$class), 
	   function(Z) names(sort(table(Z), decreasing = TRUE)), simplify = FALSE)
	   
	groupedAnalysisNew = groupedAnalysis
	pp = ncol(groupedAnalysis)
	for (j in 2:pp)
	{	
		groupedAnalysisNew[[j]] = lapply(groupedAnalysis[[j]], 
		function(Z) as.character(Z[1:(min(length(Z),maxFeatures))])) 
	}
	groupedAnalysis = groupedAnalysisNew
	 	
	groupedFrequencies = aggregate(featuresAndObs[,(2 + p.new):min((p.new + 1 + components), 2*p.new+1)], 
	list(featuresAndObs$class), function(Z) round(mean(Z),4))
	  
	colnames(groupedAnalysis)[1] = colnames(groupedFrequencies)[1] = colnames(featuresAndObs)[1] = clusterName
	newNames = vector(length = pp-1)
	for (i in 2:ncol(groupedAnalysis))	{  newNames[i-1] = paste("Component", i-1, sep = "") }
	colnames(groupedAnalysis)[-1] = newNames
	colnames(groupedFrequencies)[-1] = rm.string(colnames(groupedFrequencies)[-1], "localVariable")
	
	if (flagList)
	{
		if (!is.null(clusteredObject$classes))
		{
			featuresAndObs[,1] = clusteredObject$classes[featuresAndObs[,1]]
			groupedAnalysis[,1] = clusteredObject$classes[groupedAnalysis[,1]]
			groupedFrequencies[,1] = clusteredObject$classes[groupedFrequencies[,1]]
		}
	}
	else
	{
		if (is.factor(clusteredObject))
		{
			featuresAndObs[,1] = clusteredObject$classes[featuresAndObs[,1]]
			groupedAnalysis[,1] = clusteredObject$classes[clusterAnalysis[,1]]
			groupedFrequencies[,1] = clusteredObject$classes[groupedFrequencies[,1]]
		}
	}
	
	cat("Clustered observations:\n")
	print(head(featuresAndObs[,varIdx]))
	cat("...\n\nMost influential features by component (", maxFeatures," features per component):\n", sep="")
	print(groupedAnalysis)
	cat("\nComponent frequencies (", components, " out of ", p," possible ones):\n", sep="")
	print(groupedFrequencies)
	cat("\n\n")	
		
	if (!is.null(clusteredObject))
	{
		aggregateCategorical = aggregateNumerical = NULL
		if ((class(clusteredObject)[1]  == "unsupervised")) { Class = mergeOutliers(clusteredObject) }
		else 
		{ 
			if ( (class(clusteredObject)[1] == "randomUniformForest.formula") || (class(clusteredObject)[1] == "randomUniformForest") )
			{
				if (OOB)
				{
					if (is.null(clusteredObject$classes)) { stop("'clusteredObject' is not a classification model") }
					
					if (length(clusteredObject$forest$OOB.predicts) == n)
					{  Class = clusteredObject$classes[clusteredObject$forest$OOB.predicts] }
					else
					{  stop("OOB predictions have not their length equal to the number of rows in the data") }
				}
				else 
				{  
					if (is.null(clusteredObject$classes)) { stop("'clusteredObject' is not a classification model") }
					
					if (!is.null(clusteredObject$predictionObject))
					{	Class = clusteredObject$classes[clusteredObject$predictionObject$majority.vote]	}				
					else
					{  stop("Predictions have not their length equal to the number of rows in the data") }	
				}
			}
			else
			{ 	Class = clusteredObject  }
			Class = as.factor(Class)			
		}
		classSize = tabulate(Class)
		#names(classSize) = NULL
		if (is.null(categorical))  
		{ 
			categorical = which.is.factor(X, maxClasses = n) 
			keepIdx = which(categorical == 0)
		}
		else
		{
			if (is.character(categorical))
			{   categorical = which(colnames(X) == categorical) }
			
			keepIdx = (1:p)[-categorical]
		}
			
		pp = length(keepIdx)
		
		if (pp > 0)
		{
			localImportanceVariables = names(object$variableImportanceOverInteractions)
			keepIdx = keepIdx[rmNA(match(rmNA(match(localImportanceVariables, colnames(X))), keepIdx))]
			nFeatures = min(10, length(keepIdx))
			aggregateNumerical = aggregate(X[,keepIdx], list(Class), sum)
			aggregateNumerical = cbind(aggregateNumerical[,1], classSize, aggregateNumerical[,-1])
			colnames(aggregateNumerical)[1] = clusterName
			colnames(aggregateNumerical)[2] = "Size"
			colnames(aggregateNumerical)[-c(1,2)] = colnames(X)[keepIdx]
			cat("Numerical features aggregation (", min(10, nFeatures)," most important ones by their interactions):\n", sep="")
			if (nFeatures < 10) { print(aggregateNumerical) }
			else {  print(aggregateNumerical[,1:10]) }
			
			cat("\n")
			averageNumerical = aggregate(X[,keepIdx], list(Class), function(Z) round(mean(Z), 4))
			averageNumerical = cbind(averageNumerical[,1], classSize, averageNumerical[,-1])
			colnames(averageNumerical)[1] = clusterName
			colnames(averageNumerical)[2] = "Size"
			colnames(averageNumerical)[-c(1,2)] = colnames(X)[keepIdx]
			cat("Numerical features average (", min(10, nFeatures) ," most important ones by their interactions):\n", sep="")
			if (nFeatures < 10) { print(averageNumerical) }
			else {  print(averageNumerical[,1:10]) }
			
			cat("\n")
			standardDeviationNumerical = aggregate(X[,keepIdx], list(Class), function(Z) round(sd(Z), 4))
			standardDeviationNumerical = cbind(standardDeviationNumerical[,1], classSize, standardDeviationNumerical[,-1])
			colnames(standardDeviationNumerical)[1] = clusterName
			colnames(standardDeviationNumerical)[2] = "Size"
			colnames(standardDeviationNumerical)[-c(1,2)] = colnames(X)[keepIdx]
			cat("Numerical features (standard) deviation (", nFeatures ," most important ones by their interactions):\n", sep="")
			if (nFeatures < 10) { print(standardDeviationNumerical) }
			else {  print(standardDeviationNumerical[,1:10]) }
			
			if (pp < p)
			{
				flag = FALSE
				keepIdx2 = (1:p)[-keepIdx]
				keepIdx2 = keepIdx2[rmNA(match(rmNA(match(localImportanceVariables, colnames(X))), keepIdx2))]
				if (length(keepIdx2) == 0) {  keepIdx2 = (1:p)[-keepIdx]; flag = TRUE }
				nFeatures = min(10, length(keepIdx2))
				aggregateCategorical = aggregate(X[,keepIdx2], list(Class), function(Z) names(which.max(table(Z))))
				aggregateCategorical = cbind(aggregateCategorical[,1], classSize, aggregateCategorical[,-1])
				colnames(aggregateCategorical)[1] = clusterName
				colnames(aggregateCategorical)[2] = "Size"
				colnames(aggregateCategorical)[-c(1,2)] = colnames(X)[keepIdx2]
				if (flag)
				{	cat("\nCategorical features aggregation (", min(10, nFeatures) ," first features):\n", sep="")	}
				else
				{ 	cat("\nCategorical features aggregation (", min(10, nFeatures)," most important ones by their interactions):\n", sep="")	}
				if (nFeatures < 10) { print(aggregateCategorical) }
				else {  print(aggregateCategorical[,1:10]) }
				cat("\n")
			}
		}
			
		if (pp == 0)
		{
			keepIdx2 = (1:p)
			nFeatures = min(10, length(keepIdx2))
			aggregateCategorical = aggregate(X[,keepIdx2], list(Class), function(Z) names(which.max(table(Z))))
			aggregateCategorical = cbind(aggregateCategorical[,1], classSize, aggregateCategorical[,-1])
			colnames(aggregateCategorical)[1] = clusterName
			colnames(aggregateCategorical)[2] = "Size"
			colnames(aggregateCategorical)[-c(1,2)] = colnames(X)[keepIdx2]
			cat("Categorical features aggregation (", nFeatures ," first features):\n")
			if (nFeatures < 10) { print(aggregateCategorical) }
			else {  print(aggregateCategorical[,1:10]) }
		}
						
		# aggregateCategorical = aggregate(X[,-c(1,which(categorical != 0))], list(Class), function(Z) modX(Z))
				
		return(list(featuresAndObs = featuresAndObs, clusterAnalysis = groupedAnalysis, 
		componentAnalysis = groupedFrequencies, numericalFeaturesAnalysis = aggregateNumerical, categoricalFeaturesAnalysis = aggregateCategorical))
	}
	else
	{
		return(list(featuresAndObs = featuresAndObs, clusterAnalysis = groupedAnalysis, 
		componentAnalysis = groupedFrequencies))
	}
}


rm.coordinates <- function(object, whichOnes, seed = NULL, maxIters = NULL)
{
	Z = object$MDSModel$points = object$MDSModel$points[,-whichOnes]
	object$unsupervisedModel$centers = object$unsupervisedModel$centers[,-whichOnes] 
	
	if (is.null(seed)) { seed = as.numeric(object$params["seed"]) }
	if (is.null(maxIters)) { maxIters = 10 }
		
	if (object$params["baseModel"] == "MDShClust")
	{
		X.model <- hClust(Z, method = NULL, plotting = FALSE, k = as.numeric(object$params["clusters"]), 
		reduceClusters = FALSE, seed = seed)
		X.gapstat = NULL
		nbClusters = length(unique(X.model$cluster))
		centers = matrix(NA, nbClusters, ncol(Z))
		for (i in 1:nbClusters)
		{
			idx = which(X.model$cluster == i)
			centers[i,] = colMeans(Z[idx, ,drop = FALSE])
		}
		X.model$centers = centers
	}
	else
	{
		X.model <- kMeans(Z, NULL, k = as.numeric(object$params["clusters"]), 
		maxIters = maxIters, plotting = FALSE, algorithm = NULL, reduceClusters = FALSE, seed = seed)		
	}
	object$unsupervisedModel = X.model
	
	return(object)
}

unsupervised.randomUniformForest <- function(object, 
baseModel = c("proximity", "proximityThenDistance", "importanceThenDistance"),
endModel = c("MDSkMeans", "MDShClust", "MDS", "SpectralkMeans"),
endModelMetric = NULL,
samplingMethod = c("uniform univariate sampling", "uniform multivariate sampling", "with bootstrap"),
MDSmetric = c("metricMDS", "nonMetricMDS"),
proximityMatrix = NULL,
sparseProximities = FALSE,
outliersFilter = FALSE, 
Xtest = NULL, 
predObject = NULL, 
metricDimension = 2, 
coordinates = c(1,2),
bootstrapReplicates = 100,
clusters = NULL,
maxIters = NULL,
importanceObject = NULL,
maxInteractions = 2,
reduceClusters = FALSE, 
maxClusters = 5,
mapAndReduce = FALSE,
OOB = FALSE,
subset = NULL,
seed = 2014,
uthreads = "auto",
...)
{
	if (!is.null(predObject))
	{
		if (is.null(predObject$votes.data))
		{ stop("predObject must be provided using option type = 'all' when calling the predict() function.\n") }
	}	
	
	if (!is.null(Xtest)) 
	{ 	
		Xtest = NAfactor2matrix(Xtest, toGrep = "anythingParticular")
		if (length(which(is.na(Xtest)) > 0) ){ stop("NA found in data. Please treat them outside of the model.\n") }
	}
	
	flagBig = FALSE
	objectClass = class(object)
	if ( all(objectClass != "randomUniformForest") )
	{
		if (!is.null(subset)) { object = object[subset,] }		
		X = object
		n = nrow(X)
		
		if ( (n > 10000) | mapAndReduce) 
		{
			set.seed(seed)
			subsetIdx = sample(n, min(10000, floor(n/2)))
			bigX = X[-subsetIdx, ]
			X = X[subsetIdx,]
			n = nrow(X)
			flagBig  = TRUE
		}
		
		if (samplingMethod[1] == "with bootstrap")
		{ 
			cat("'with bootstrap' is only needed as the second argument of 'method' option. Default option will be computed.\n")
			samplingMethod[1] = "uniform univariate sampling"
		}
		XY <- unsupervised2supervised(X, method = samplingMethod[1], seed = seed, 
			bootstrap = if (length(samplingMethod) > 1) { TRUE } else {FALSE })
		
		cat("Created synthetic data giving them labels.\n") 
		if ((n > 10000) | mapAndReduce)
		{	
			rUF.model <- rUniformForest.big(XY$X, as.factor(XY$Y), BreimanBounds = FALSE, unsupervised = TRUE, unsupervisedMethod = samplingMethod[1], randomCut = TRUE,...)  
		}
		else
		{ 	
			rUF.model <- randomUniformForest(XY$X, as.factor(XY$Y), BreimanBounds = FALSE, 
			unsupervised = TRUE, unsupervisedMethod = samplingMethod[1],...)  
		}
		
		if (is.null(Xtest)) { Xtest = XY$X[1:n,] }
	}
	else
	{	
		if (any(objectClass == "randomUniformForest")) { rUF.model = object	}
		else { stop("Data are missing or object is not of class randomUniformForest.") }
		n = if (is.null(object$unsupervised)) { length(object$y) } else {  length(object$y)/2 }
	}
		
	if (!is.null(Xtest)) 
	{ 
		n = nrow(Xtest) 
		if ( ((n > 10000) | mapAndReduce) & !flagBig ) 
		{
			set.seed(seed)
			subsetIdx = sample(n, min(10000, floor(n/2)))
			bigX = Xtest[-subsetIdx, ]
			Xtest = Xtest[subsetIdx,]
			n = nrow(Xtest)
			flagBig  = TRUE
		}		
	}	
			
	if (baseModel[1] != "importanceThenDistance")
	{
		if (is.null(proximityMatrix))
		{
			proxMat <- proximitiesMatrix(rUF.model, fullMatrix = TRUE, Xtest = Xtest, predObject = predObject, 
			sparseProximities = sparseProximities, pthreads = uthreads) 
		}
		else
		{ proxMat = proximityMatrix }
	}
	else
	{
		if (!is.null(importanceObject)) { imp.rUFModel = importanceObject }
		else { imp.rUFModel <- importance(rUF.model, Xtest = Xtest, maxInteractions = maxInteractions) }
		proxMat <- observationsImportance(Xtest, imp.rUFModel) 
	}
	
	p = ncol(proxMat)
	
	grepMDS = length(grep("MDS", endModel[1]))
	grepSpectral = length(grep("Spectral", endModel[1]))
	if ((grepMDS == 1) | (grepSpectral == 1))
	{
		if (grepSpectral == 1)
		{
			X.scale = list()
			cat("Spectral clustered points have been put on ...$MDSModel for compatibility with new points and others clustering objects that have to be updated.\n\n")
			Z <- specClust(proxMat, k = max(3, metricDimension))
			Z = Z[, coordinates]			
		}
		else
		{
			distanceProxy = TRUE
			if (baseModel[1] == "proximity") { 	distanceProxy = FALSE }		
			X.scale <- MDSscale(proxMat, metric = MDSmetric[1], dimension = metricDimension, distance = distanceProxy, plotting = FALSE, seed = seed)
			
			Z = X.scale$points[, ,drop = FALSE]
			
			if (dim(Z)[2] == 0) 
			stop ("MDS can not be achieved. Data can probably be not clustered using this (randomUniformForest) dissimilarity matrix.\nOptions used might be reconsidered.\n")
		}

		if (flagBig)
		{
			pZ = ncol(Z)
			X.scale.ruf = vector('list', pZ)
			newZ = matrix(NA, nrow(Z) + nrow(bigX), pZ)
			if (mapAndReduce)
			{
				cat("Entered in regression mode to learn MDS/spectral points.\n")
				for (i in 1:pZ)	
				{ 
					X.scale.ruf[[i]] = rUniformForest.big(X, Z[,i], randomCut = TRUE,...) 	
					newZ[-subsetIdx,i] = predict(X.scale.ruf[[i]], bigX)
					newZ[subsetIdx,i] = Z[,i]
				}
			}
			else
			{
				for (i in 1:pZ)	
				{ 
					X.scale.ruf[[i]] = randomUniformForest(X, Z[,i], ...) 
					newZ[-subsetIdx,i] = predict(X.scale.ruf[[i]], bigX)
					newZ[subsetIdx,i] = Z[,i]
				}			
			}
			Z = newZ
		}
				
		nn = nrow(Z)
		if (endModel[1] == "MDS") { outliersFilter = FALSE }
		if (outliersFilter)
		{
			idx1 = which( (Z[,1] < quantile(Z[,1], 0.025)) & (Z[,2] < quantile(Z[,2], 0.025)) ) 
			idx2 = which( (Z[,1] > quantile(Z[,1], 0.975)) & (Z[,2] > quantile(Z[,2], 0.975)))
			idx12 = c(idx1, idx2)
			if (length(idx12) > 0)	
			{ 
				cat("Outliers have been removed for the two first coordinates.\n")
				Z = Z[-idx12,]  
			}
			nn = nrow(Z)
		}
		
		if (is.null(predObject))
		{
			if (endModel[1] == "MDS")
			{
				X.gapstat = NULL 
				if (is.null(object$unsupervised))
				{
					if (OOB)
					{
						if (!is.null(object$forest$OOB.predicts))
						{	predictions = object$forest$OOB.predicts }
						else
						{ 
							cat("OOB option has been called, but there is currently no OOB classifier. Hence true labels have been used.\n") 
							predictions = object$y
						}						
					}
					else
					{	
						if (is.null(object$predictionObject)) { predictions = object$y }
						else { 	predictions = object$predictionObject$majority.vote }				
					}
					
					classes = sort(unique(predictions))
					k = length(classes)
					
					{  
						centers = matrix(0, k, ncol(Z))
						for (i in 1:k)
						{	
							idx = which(predictions == classes[i])
							centers[i,] = colMeans(Z[idx,])
						}						
					}
					
					rownames(centers) = as.character(classes)
					X.model = list(cluster = predictions, centers = centers)
				}
				else
				{
					cat("Clustering can not be done except in the supervised mode.\nPlease send ...$MDSModel$points to a clustering algorithm.\n")
					X.model = list(cluster = NULL, centers = NULL)
				}
			}
					
			if ( (endModel[1] == "MDSkMeans") | (endModel[1] == "SpectralkMeans") )
			{
				if (is.null(clusters))  
				{	X.gapstat <- gap.stats(Z, B = bootstrapReplicates, maxClusters = maxClusters, seed = seed) 	}
				else 
				{ X.gapstat = NULL }
				
				X.model <- kMeans(Z, X.gapstat, k = clusters, maxIters = maxIters, plotting = FALSE, algorithm = endModelMetric, reduceClusters = reduceClusters,  seed = seed)

				if (endModel[1] == "SpectralkMeans") { X.scale$points = Z }
			}
			
			if (endModel[1] == "MDShClust")  
			{	
				X.model <- hClust(Z, method = endModelMetric, plotting = FALSE, k = clusters, 
				reduceClusters = reduceClusters, seed = seed)
				X.gapstat = NULL
				nbClusters = length(unique(X.model$cluster))
				centers = matrix(NA, nbClusters, ncol(Z))
				for (i in 1:nbClusters)
				{
					idx = which(X.model$cluster == i)
					centers[i,] = colMeans(Z[idx,, drop = FALSE])
				}
				X.model$centers = centers
			}

			if (outliersFilter)
			{
				if (length(idx12) > 0)
				{
					if ( (endModel[1] == "MDSkMeans") | (endModel[1] == "SpectralkMeans") ) 
					{	clusterOutliers <- which.is.nearestCenter(Z[c(idx1,idx2),], X.model$centers)	}
					if (endModel[1] == "MDShClust")  
					{	
						nbClusters = length(unique(X.model$cluster))
						centers = matrix(NA, nbClusters, ncol(Z))
						for (i in 1:nbClusters)
						{
							idx = which(X.model$cluster == i)
							centers[i,] = colMeans(Z[idx,, drop = FALSE])
						}
						
						clusterOutliers <- which.is.nearestCenter(Z[c(idx1,idx2),], centers)
					}
					names(clusterOutliers) = c(idx1,idx2)
					X.model$clusterOutliers = clusterOutliers
				}
			}
		}
		else
		{
			predictions = predObject$majority.vote
			classes = sort(unique(predictions))
			k = length(classes)
			centers = matrix(0, k, ncol(Z))
			for (i in 1:p)
			{	
				idx = which(predictions == classes[i])
				centers[i,] = colMeans(Z[idx,])
			}
			rownames(centers) = as.character(classes)
			X.model = list(cluster = predictions, centers = centers)
			X.gapstat = NULL 
		}
	}
	else
	{
		if (endModel[1] == "kMeans")
		{
			X.scale = NULL
			if (is.null(clusters))  
			{	
				X.gapstat <- gap.stats(proxMat[sample(n, floor(n/2)),, drop = FALSE], B = bootstrapReplicates, maxClusters = maxClusters,  seed = seed) 
			}
			else 
			{ X.gapstat = NULL }		
			X.model <- kMeans(proxMat, X.gapstat, k = clusters, maxIters = maxIters, plotting = FALSE, algorithm = endModelMetric,  seed = seed)
		}
		
		if (endModel[1] == "hClust")  
		{	
			X.scale = NULL
			X.model <- hClust(proxMat, method = endModelMetric, plotting = FALSE, k = clusters,  seed = seed)
			X.gapstat = NULL 
		}
	}
	clusters = length(unique(X.model$cluster))	
		
	modelParams = c(baseModel[1], endModel[1],  if (is.null(endModelMetric)) { "NULL" } else { endModelMetric }, 
	samplingMethod[1], MDSmetric[1], is.null(Xtest), is.null(predObject), metricDimension, concatCore(as.character(coordinates)), bootstrapReplicates, clusters, if (is.null(maxIters)) { "NULL" } else { maxIters }, 
	maxInteractions, maxClusters, mapAndReduce, outliersFilter, reduceClusters, seed, sparseProximities)

	names(modelParams) = c("baseModel", "endModel", "endModelMetric", "samplingMethod", "MDSmetric", "Xtest", "predObject", "metricDimension", "coordinates", "bootstrapReplicates", "clusters", "maxIters", "maxInteractions", "maxClusters", "mapAndReduce", "outliersFilter", "reduceClusters", "seed", "sparseProximities")
		
	if (flagBig) 
	{ 	
		X.scale$points = Z
		largeDataLearningModel = X.scale.ruf
		unsupervisedObject = list(proximityMatrix = proxMat, MDSModel = X.scale, unsupervisedModel = X.model, largeDataLearningModel = largeDataLearningModel, gapStatistics = X.gapstat, rUFObject = rUF.model, 
		nbClusters = clusters, params = modelParams) 	
	}
	else
	{
		unsupervisedObject = list(proximityMatrix = proxMat, MDSModel = X.scale, unsupervisedModel = X.model, gapStatistics = X.gapstat, rUFObject = rUF.model, nbClusters = clusters, params = modelParams) 	
	}
	
	class(unsupervisedObject) <- "unsupervised"
	
	unsupervisedObject
}
# END OF FILE