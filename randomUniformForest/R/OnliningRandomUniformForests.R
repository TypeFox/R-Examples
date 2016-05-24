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

onlineCombineRUF <- function(rufObject1, rufObject2)
{
	rufObject.0 <- mergeLists(rufObject1, rufObject2, OOBEval = TRUE)	
	rufObject.0$regression = rufObject1$regression	
	if (!is.null(rufObject.0$percent.varExplained))
	{	rufObject.0$percent.varExplained = c(rufObject1$percent.varExplained, rufObject2$percent.varExplained)	}	
	if (!is.null(rufObject.0$OOB))
	{	
		if (is.character(rufObject.0$OOB))
		{	rufObject.0$OOB = rufObject1$OOB }
		else
		{
			p1 = ncol(rufObject1$OOB)
			p2 = ncol(rufObject2$OOB)
			tmp.OOB1 = rufObject1$OOB[,-p1]
			tmp.OOB2 = rufObject2$OOB[,-p2]
						
			if (p1 != p2)
			{
				if (p1 > p2)
				{
					tmp.OOB = tmp.OOB1
					matchIdx = rmNA(match(as.numeric(colnames(tmp.OOB2)), as.numeric(colnames(tmp.OOB1))))
					for (i in 1:length(matchIdx))
					{	
						tmp.OOB[matchIdx[i],matchIdx[i]] = tmp.OOB[matchIdx[i],matchIdx[i]] + tmp.OOB2[matchIdx[i],matchIdx[i]]
					}	
				}
				else
				{
					tmp.OOB = tmp.OOB2
					matchIdx = rmNA(match(as.numeric(colnames(tmp.OOB1)), as.numeric(colnames(tmp.OOB2))))
					for (i in 1:length(matchIdx))
					{	
						tmp.OOB[matchIdx[i],matchIdx[i]] = tmp.OOB[matchIdx[i],matchIdx[i]] + tmp.OOB1[matchIdx[i],matchIdx[i]]
					}				
				}
				class.error =  1 - diag(tmp.OOB)/rowSums(tmp.OOB)	
				rufObject.0$OOB = cbind(tmp.OOB, class.error)
			}
			else
			{			
				tmp.OOB = tmp.OOB1 + tmp.OOB2
				class.error =  1 - diag(tmp.OOB)/rowSums(tmp.OOB)	
				rufObject.0$OOB = cbind(tmp.OOB, class.error)
			}
		}
	}	
	tmp.rufObject.0 = rufObject.0	
	if (!is.null(rufObject.0$OOB.strengthCorr))
	{	
		rufObject.0$OOB.strengthCorr = NULL		
		if (!is.null(rufObject1$OOB.strengthCorr$PE.forest))
		{
			rufObject.0$OOB.strengthCorr$PE.forest = c(rufObject1$OOB.strengthCorr$PE.forest, rufObject2$OOB.strengthCorr$PE.forest)
			rufObject.0$OOB.strengthCorr$PE.max = c(rufObject1$OOB.strengthCorr$PE.max, rufObject2$OOB.strengthCorr$PE.max)
			rufObject.0$OOB.strengthCorr$PE.tree = c(rufObject1$OOB.strengthCorr$PE.tree, rufObject2$OOB.strengthCorr$PE.tree)
			rufObject.0$OOB.strengthCorr$mean.corr = c(rufObject1$OOB.strengthCorr$mean.corr, rufObject2$OOB.strengthCorr$mean.corr)
		}
		else
		{
			rufObject.0$OOB.strengthCorr$PE = c(rufObject1$OOB.strengthCorr$PE, rufObject2$OOB.strengthCorr$PE)
			rufObject.0$OOB.strengthCorr$avg.corr = c(rufObject1$OOB.strengthCorr$avg.corr, rufObject2$OOB.strengthCorr$avg.corr)
			rufObject.0$OOB.strengthCorr$strength = c(rufObject1$OOB.strengthCorr$strength, rufObject2$OOB.strengthCorr$strength)
			rufObject.0$OOB.strengthCorr$std.strength = c(rufObject1$OOB.strengthCorr$std.strength, rufObject2$OOB.strengthCorr$std.strength)
		}
	}	
	if (!is.null(rufObject.0$variableImportance))
	{	
		rufObject.0$variableImportance  = NULL
		tmp.var.imp.names = unique(c(as.character(rufObject1$variableImportance[,1]), as.character(rufObject2$variableImportance[,1])))
		p = length(tmp.var.imp.names)
		score = classFrequency = classEstimate = rep(0,p)		
		for (i in 1:p)
		{
			idx.1 = which(rufObject1$variableImportance[,1] == tmp.var.imp.names[i])
			idx.2 = which(rufObject2$variableImportance[,1] == tmp.var.imp.names[i])		    
			if (length(idx.1) > 0)
			{	
				score[i] = rufObject1$variableImportance[idx.1,2]				
				if (!rufObject.0$regression) 
				{
					if (any(is.na(rufObject1$variableImportance[,3])))
					{	rufObject1$variableImportance[,3] = as.vector(na.replace(matrix(rufObject1$variableImportance[,3]))) }
					classFrequency[i] = rufObject1$variableImportance[idx.1,4]
					classEstimate[i] = rufObject1$variableImportance[idx.1,3]
				}
			}			
			if  (length(idx.2) > 0)
			{	
				score[i] = score[i]  + rufObject2$variableImportance[idx.2, 2]				
				if (!rufObject.0$regression)
				{
					if (any(is.na(rufObject2$variableImportance[,3])))
					{	rufObject2$variableImportance[,3] = as.vector(na.replace(matrix(rufObject2$variableImportance[,3]))) }
					if (classEstimate[i] == rufObject2$variableImportance[idx.2,3])
					{	classFrequency[i] = 0.5*(classFrequency[i] + rufObject2$variableImportance[idx.2,4]) }
					else
					{
						if (!is.na(classFrequency[i]))	
						{	allClassFrequency = c(classFrequency[i], rufObject2$variableImportance[idx.2,4])	}
						else
						{	 allClassFrequency = rufObject2$variableImportance[idx.2,4]	 }						
						classFrequency[i] = allClassFrequency[which.max(allClassFrequency)]
						classEstimate[i] = if( allClassFrequency[1] == 0) { rufObject1$variableImportance[idx.2,3] }
						else { c(classEstimate[i], rufObject1$variableImportance[idx.2,3])[which.max(allClassFrequency)] } 
					}
				}
			}  
		}		
		if (!rufObject.0$regression)
		{
			var.imp.output = cbind(score, classEstimate, classFrequency, round(100*score/max(score),2), round(100*score/sum(score),0))
			var.imp.output =  data.frame(tmp.var.imp.names, var.imp.output)
			var.imp.output = var.imp.output[order(var.imp.output[,2], decreasing = TRUE),]					
			colnames(var.imp.output) = c("variables", "score", "class", "class.frequency", "percent", "percent.importance")
		}
		else
		{
			var.imp.output = cbind(score,  round(100*score/max(score),2), round(100*score/sum(score),0))
			var.imp.output =  data.frame(tmp.var.imp.names, var.imp.output)
			var.imp.output = var.imp.output[order(var.imp.output[,2], decreasing = TRUE),]					
			colnames(var.imp.output) = c("variables", "score", "percent", "percent.importance")
		}		
		rufObject.0$variableImportance = var.imp.output
	}
	if (!is.null(rufObject.0$OOB.strengthCorr)) { 	tmp.rufObject.0$OOB.strengthCorr = rufObject.0$OOB.strengthCorr	}	
	if (!is.null(rufObject.0$variableImportance))	{   tmp.rufObject.0$variableImportance  = rufObject.0$variableImportance  }	
	return(tmp.rufObject.0)
}

rUniformForest.big <- function(X, Y = NULL, xtest = NULL, ytest = NULL, nforest = 2, randomCut = FALSE, 
	reduceDimension = FALSE,
    reduceAll = FALSE, 			
	replacement = FALSE, 
	subsample = FALSE,	
	ntree = 100,
    nodesize = 1, 
	maxnodes = Inf,
	mtry = ifelse(bagging,ncol(X),floor(4/3*ncol(X))),
	regression = ifelse(is.factor(Y),FALSE, TRUE),
    subsamplerate = ifelse(regression, 0.7, 1),
	replace = ifelse(regression,FALSE,TRUE),
	OOB = TRUE,
	BreimanBounds = ifelse(OOB, TRUE, FALSE),
    depth = Inf,
    depthcontrol = NULL,
    importance = TRUE,	
	bagging = FALSE,
	unsupervised = FALSE, 
	unsupervisedMethod = c("uniform univariate sampling", "uniform multivariate sampling", "with bootstrap"),
	classwt = NULL,	
	oversampling = 0,
	targetclass = -1,
	outputperturbationsampling = FALSE,
    rebalancedsampling = FALSE,
	featureselectionrule = c("entropy", "gini", "random", "L1", "L2"),	
	randomcombination = 0,
	randomfeature = FALSE,
	categoricalvariablesidx = NULL,
	na.action = c("fastImpute", "accurateImpute", "omit"),
	logX = FALSE,
	classcutoff = c(0,0),
	subset = NULL,
	usesubtrees = FALSE,
	threads = "auto",
	parallelpackage = "doParallel")
{	
	if (!is.null(subset))
	{
		X = X[subset,]
		Y = Y[subset]
	}
	
	if (exists("categorical", inherits = FALSE)) {  categoricalvariablesidx = categorical }
	else { categorical = NULL  }
	if ((nrow(X)/nforest < 50) & !reduceDimension & !reduceAll) 
	{ stop("Too small size to achieve efficiency in learning")	}
	X <- fillVariablesNames(X)
	if (is.null(Y)) 
	{ 
		cat("Enter in unsupervised learning.\n")
		if (unsupervisedMethod[1] == "with bootstrap")
		{ 
			cat("'with bootstrap' is only needed as the second argument of 'method' option. Default option will be computed.\n")
			unsupervisedMethod[1] = "uniform univariate sampling"
		}
		XY <- unsupervised2supervised(X, method = unsupervisedMethod[1], bootstrap = if (length(unsupervisedMethod) > 1) { TRUE } else {FALSE })
		X = XY$X
		Y = as.factor(XY$Y)
		unsupervised = TRUE
	}	
    getFactors <- which.is.factor(X, count = TRUE) 	
	if (is.data.frame(X)) {  cat("X is a dataframe. String or factors have been converted to numeric values.\n")  }
	X <- NAfactor2matrix(X)	
	if (randomCut)
	{
		set.seed(sample(1e9,1))
		n = nrow(X)
		newIdx = sample(n,n)
		X = X[newIdx,]
		Y = Y[newIdx]
	}
	if ( length(which( (X == Inf) | (X == -Inf) ) ) > 0) 
	{ stop("Inf or -Inf values found in data. Learning can not be done.\nPlease replace them with NA in order to learn") }
	XY <- NATreatment(X, Y, na.action = na.action, regression = regression)
	X  <- XY$X
	Y  <- XY$Y
	rm(XY)
	if (randomcombination[1] > 0)
	{ 
		randomcombinationString = randomcombination[1]
		L.combination = length(randomcombination)
		if (L.combination%%3 != 0)
		{	weights = round(replicate(length(randomcombination)/2, runif(1)), 2)	}
		else
		{ 	
			weights = round(randomcombination[(L.combination + 1 - L.combination/3):L.combination],2)
			randomcombination = randomcombination[1:(L.combination - (L.combination/3))]
		}		
		for (i in 2:length(randomcombination)) {  randomcombinationString = paste(randomcombinationString, randomcombination[i], sep=",") }
		for (i in 1:length(weights)) {  randomcombinationString = paste(randomcombinationString, weights[i], sep=",") }
		if (!is.null(xtest)) { xtest <- randomCombination(xtest, combination = randomcombination, weights = weights) }
		randomcombinationObject = list(randomcombination, weights)
		X <- randomCombination(X, combination = randomcombinationObject[[1]], weights = randomcombinationObject[[2]])	
	}
	else
	{   randomcombinationObject = randomcombination }
	RUF.model <- randomUniformForestCore.big(X, Y, nforest, reduceDimension = reduceDimension, reduceAll = reduceAll, r.ntree = ntree, r.nodeMinSize = nodesize, r.maxNodes = maxnodes, r.depth = depth, r.depthControl = depthcontrol, r.features = mtry, r.bootstrap = replace, r.treeSubsampleRate = subsamplerate, r.overSampling = oversampling, r.outputPerturbation = outputperturbationsampling, r.targetClass = targetclass, 
	r.rebalancedSampling = rebalancedsampling, r.classcutoff = classcutoff, r.use.OOB = OOB,
	r.BreimanBounds = BreimanBounds, r.treeBagging = bagging, r.randomCombination = randomcombination, 
	r.randomFeature = randomfeature, r.variableImportance = importance, r.whichCatVariables = categoricalvariablesidx, r.featureSelectionRule = featureselectionrule, r.logX = logX, r.classwt = classwt, r.regression = regression,  r.useSubTrees = usesubtrees, r.threads = threads, r.parallelPackage = parallelpackage, replacement = FALSE, subsample = FALSE)		
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
		{ 	rebalancedsamplingString = paste(rebalancedsamplingString, rebalancedsampling[i], sep=",") }	
	}	
	if (RUF.model$regression)	{  classcutoff = c(0,0)  }	
	if (as.numeric(classcutoff[2]) == 0)	{  classcutoffString = FALSE }
	else
	{
		classcutoffString = levels(Y)[which(levels(Y) == as.character(classcutoff[1]))]
		classcutoffString = paste("Class ", classcutoffString, "," , as.numeric(classcutoff[2])*100, "%", sep ="")
	}	
	paramsObject = c(reduceDimension, reduceAll, nforest, ntree, mtry, nodesize, maxnodes, as.character(replace), 		bagging, depth, 
		ifelse(is.null(depthcontrol), length(depthcontrol)> 0, depthcontrol), OOB, importance, subsamplerate, 
		ifelse(is.null(classwt),  length(classwt) > 0, classwtString), classcutoffString, 
		ifelse((oversampling == 0),(oversampling != 0), oversampling), outputperturbationsampling, 
		ifelse(length(targetclass) > 1, targetclassString, targetclass), 
		ifelse(length(rebalancedsampling) > 1, rebalancedsamplingString, rebalancedsampling), 
		ifelse((randomcombination[1] == 0),(randomcombination != 0), randomcombinationString), randomfeature, 
		ifelse(is.null(categoricalvariablesidx), length(categoricalvariablesidx) > 0, length(categoricalvariablesidx)), featureselectionrule[1])	
	names(paramsObject) = c("reduceDimension", "reduceAll", "nforest", "ntree", "mtry", "nodesize", "maxnodes", "replace", "bagging", "depth", "depthcontrol",  "OOB", "importance", "subsamplerate", "classwt", "classcutoff", "oversampling", "outputperturbationsampling", "targetclass", "rebalancedsampling", "randomcombination", "randomfeature", "categorical variables", "featureselectionrule")	
	if (RUF.model$regression) 
	{ 
		paramsObject[24] = "Sum of squared residuals"
	    if (paramsObject[8] == "FALSE") { paramsObject[11] = FALSE }
	}		
	paramsObject = as.data.frame(paramsObject)
	RUF.model$logX = logX	
	RUFObject <- genericOutput(xtest, ytest, paramsObject, RUF.model, ytrain = Y, classcutoff = classcutoff)	
	RUFObject$logX = logX
	if(is.null(colnames(X))) 
	{  
		varNames = NULL
		for (i in 1:ncol(X)) { varNames = c(varNames,paste("V", i, sep="")) }
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
	if (!is.null(Y)) { 	RUFObject$y = Y	 }
	if (!is.null(categoricalvariablesidx ))
	{  RUFObject$categoricalvariables = categoricalvariablesidx }
	if (unsupervised) 
	{ 
		RUFObject$unsupervised = TRUE 
		RUFObject$unsupervisedMethod  = if (length(unsupervisedMethod) == 3) 
										{ unsupervisedMethod[1] } 
										else 
										{ 
											if (length(unsupervisedMethod) == 1) { unsupervisedMethod } 
											else { unsupervisedMethod[1:2] } 
										}
	}
	RUFObject$call <- match.call()
	class(RUFObject) <- "randomUniformForest"	
	RUFObject
}

randomUniformForestCore.big <- function(X, Y, howManyForests, reduceDimension = FALSE, reduceAll = FALSE, 
	replacement = FALSE, 
	subsample = FALSE,
	r.nodeMinSize = 1,
	r.maxNodes = Inf, 
	r.depth = Inf,
	r.depthControl = NULL,
	r.features = floor(4/3*ncol(X)), 
	r.bootstrap = ifelse(r.regression, FALSE, TRUE), 
	r.overSampling = 0,
	r.outputPerturbation = FALSE,
	r.targetClass = 0,
	r.rebalancedSampling = FALSE,
	r.regression = TRUE, 
	r.ntree = 100, 
	r.use.OOB = TRUE,
	r.BreimanBounds = ifelse(r.use.OOB, TRUE, FALSE),
	r.treeSubsampleRate = ifelse(r.regression, 0.7, 1),
	r.variableImportance = TRUE,
	r.treeBagging = FALSE,
	r.classwt = NULL,
	r.classcutoff = c(0,0),
	r.randomCombination = 0,
	r.randomFeature = FALSE,
	r.whichCatVariables = NULL,
	r.featureSelectionRule = c("entropy", "gini", "random", "L1", "L2"),
	r.logX = FALSE,
	r.useSubTrees = FALSE,
	r.threads = "auto",
	r.parallelPackage = "doParallel")
{	
	nCovForests = -1
	r.ntree = round(r.ntree/(howManyForests + nCovForests + 1),0)		
	if (howManyForests > 1)
	{
		if (reduceAll)
		{
			X.list_kdimensions <- setManyDatasets(X, howManyForests, dimension = TRUE, replace = replacement, 
			subsample = subsample)
			X.list_kdatasets <- setManyDatasets(X, howManyForests, replace = replacement, subsample = subsample)		
			X.set = X[X.list_kdatasets[[1]], X.list_kdimensions[[1]]]
			Y.set = Y[X.list_kdatasets[[1]]] 
		}
		else
		{
			if (reduceDimension)
			{
				X.list_kdimensions <- setManyDatasets(X, howManyForests, dimension = reduceDimension, replace = replacement, subsample = subsample)
				X.set = X
				X.set[,-X.list_kdimensions[[1]]] = 0
				Y.set = Y
			}
			else
			{
				X.list_kdatasets <- setManyDatasets(X, howManyForests, replace = replacement, subsample = subsample)
				X.set = X[X.list_kdatasets[[1]],]
				Y.set = Y[X.list_kdatasets[[1]]] 
			}
		}
	}
	else
	{ X.set = X;  Y.set = Y }
	
	rUF1 <- randomUniformForestCore(X.set, trainLabels = Y.set, depth = r.depth, depthControl = r.depthControl, 
		ntree = r.ntree, nodeMinSize = r.nodeMinSize, maxNodes = r.maxNodes, features = r.features, 
		rf.bootstrap = r.bootstrap, use.OOB = r.use.OOB, BreimanBounds = r.BreimanBounds, threads = r.threads, rf.treeSubsampleRate = r.treeSubsampleRate, rf.overSampling = r.overSampling, rf.targetClass = r.targetClass, rf.rebalancedSampling = r.rebalancedSampling, variableImportance = r.variableImportance, 
		rf.treeBagging = r.treeBagging, rf.regression = r.regression, rf.randomCombination = r.randomCombination, 
		rf.randomFeature = r.randomFeature, whichCatVariables = r.whichCatVariables, logX = r.logX, classwt = r.classwt, 
		classCutOff = r.classcutoff, rf.outputPerturbationSampling = r.outputPerturbation, rf.featureSelectionRule = r.featureSelectionRule, useSubTrees = r.useSubTrees, parallelPackage = r.parallelPackage[1])
	
	if (howManyForests > 1)
	{
		for (i in 2:(howManyForests + nCovForests + 1))
		{
			{
				if (reduceAll)
				{	
					X.set = X[X.list_kdatasets[[i]],X.list_kdimensions[[i]]]
					Y.set = Y[X.list_kdatasets[[i]]] 
				}
				else
				{
					if (reduceDimension) 
					{ 
						X.set = X
						X.set[,-X.list_kdimensions[[i]]] = 0
					}
					else
					{
						X.set = X[X.list_kdatasets[[i]],]
						Y.set = Y[X.list_kdatasets[[i]]] 
					}
				}
			}
			rUF2 <- randomUniformForestCore(X.set, trainLabels = Y.set, depth = r.depth, depthControl = r.depthControl, ntree = r.ntree, nodeMinSize = r.nodeMinSize, maxNodes = r.maxNodes, features = r.features, 
			rf.bootstrap = r.bootstrap, use.OOB = r.use.OOB, BreimanBounds = r.BreimanBounds, threads = r.threads, rf.treeSubsampleRate = r.treeSubsampleRate, rf.overSampling = r.overSampling, 
			rf.targetClass = r.targetClass, rf.rebalancedSampling = r.rebalancedSampling,
			variableImportance = r.variableImportance, rf.treeBagging = r.treeBagging, rf.regression = r.regression, rf.randomCombination = r.randomCombination, rf.randomFeature = r.randomFeature, whichCatVariables = r.whichCatVariables, logX = r.logX, classwt = r.classwt, classCutOff = r.classcutoff, rf.outputPerturbationSampling = r.outputPerturbation, rf.featureSelectionRule = r.featureSelectionRule, useSubTrees = r.useSubTrees, parallelPackage = r.parallelPackage[1])			
			rUF1 <- onlineCombineRUF(rUF1, rUF2)
		}
	}	
	return(rUF1)
}

rUniformForest.merge <- function(X = NULL, Y = NULL, xtest = NULL, ytest = NULL, previousObject = NULL, 
	newObject = NULL,
	ntree = 100, 
	mtry = ifelse(bagging,ncol(X),floor(4/3*ncol(X))),
	nodesize = 1,
	maxnodes = Inf, 
    subsamplerate = 1,
    depth = Inf,
    depthcontrol = NULL,	    
	regression = ifelse(is.factor(Y),FALSE,TRUE),
	replace = ifelse(regression,FALSE,TRUE),
	OOB = ifelse(replace,TRUE,FALSE),
	BreimanBounds = ifelse(OOB, TRUE, FALSE),
    importance = TRUE,	
	bagging = FALSE,
	classwt = NULL,	
	unsupervised = FALSE, 
	unsupervisedMethod = c("uniform univariate sampling", "uniform multivariate sampling", "with bootstrap"),	
	oversampling = 0,
	targetclass = -1,
	outputperturbationsampling = FALSE,
    rebalancedsampling = FALSE,
	featureselectionrule = c("entropy", "gini", "random", "L1", "L2"),	
	randomcombination = 0,
	randomfeature = FALSE,
	categoricalvariablesidx = NULL,
	threads = "auto",
	logX = FALSE,
	classcutoff = c(0,0),
	subset = NULL,
	usesubtrees = FALSE,
	na.action = c("fastImpute", "accurateImpute", "omit"),
	parallelpackage = "doParallel"  )
{
	if (exists("categorical", inherits = FALSE)) {  categoricalvariablesidx = categorical }
	else { categorical = NULL  }
	
	if (!is.null(X) & !is.null(Y))
	{
		if (!is.null(subset))
		{
			X = X[subset,]
			Y = Y[subset]
		}		
		if (ntree < 2) { stop("Please use at least 2 trees for computing forest\n") }
		if ( (subsamplerate == 1) &  (replace == FALSE))  { OOB = FALSE }
		if (depth < 3) { stop("Stumps are not allowed. Minimal depth is 3, leading to 4 leaf nodes.\n") }
		if (BreimanBounds &  ntree > 500) { cat("Warning : Breiman Bounds, especially for multiclass problems, are computationnaly intensive.\n") }
		if (maxnodes < 6) { stop("Maximal number of nodes must be above 5.\n") }
		{
			X <- fillVariablesNames(X)			
			getFactors <- which.is.factor(X, count = TRUE)		 
			if (is.data.frame(X)) 
			{ 
				cat("X is a dataframe. String or factors have been converted to numeric values.\n") 
				X <- NAfactor2matrix(X)	
			}					
			if (!is.null(categoricalvariablesidx))
			{
				if (categoricalvariablesidx[1] == "all")
				{	
					factorVariables <- which(getFactors > 0)
					if (length(factorVariables) > 0) 	{ 	categoricalvariablesidx = factorVariables	}
					else
					{ 
					  cat("\nNo categorical variables found. Please type them manually, replacing in categoricalvariablesidx
					  option argument 'all' by a vector of the variables that need to be treated as categorical.\n")
					}
				}
				else 
				{  
					if (is.character(categoricalvariablesidx[1]))
					{	categoricalvariablesidx = sort(match(categoricalvariablesidx, colnames(X)))	}
				}
			}
			if ( length(which( (X == Inf) | (X == -Inf) ) ) > 0) 
			{ stop("Inf or -Inf values found in data. Learning can not be done.\nPlease replace them with NA in order to learn") }
			XY <- NATreatment(X, Y, na.action = na.action, regression = regression)
			X  <- XY$X
			Y  <- XY$Y
			rm(XY)
			if (randomcombination[1] > 0)
			{ 
				randomcombinationString = randomcombination[1]
				L.combination = length(randomcombination)
				if (L.combination%%3 != 0)
				{	weights = round(replicate(length(randomcombination)/2, runif(1)), 2)	}
				else
				{ 	
					weights = round(randomcombination[(L.combination + 1 - L.combination/3):L.combination],2)
					randomcombination = randomcombination[1:(L.combination - (L.combination/3))]
				}				
				for (i in 2:length(randomcombination)) 
				{  randomcombinationString = paste(randomcombinationString, randomcombination[i], sep=",") }			
				for (i in 1:length(weights)) 
				{  randomcombinationString = paste(randomcombinationString, weights[i], sep=",") }				
				if (!is.null(xtest)) { xtest <- randomCombination(xtest, combination = randomcombination, weights = weights) }				
				randomcombinationObject = list(randomcombination, weights)						
				X <- randomCombination(X, combination = randomcombinationObject[[1]], 
				weights = randomcombinationObject[[2]])				
			}
			else
			{ randomcombinationObject = randomcombination }
		}
	}
	RUF.model <- randomUniformForestCore.merge(X = X, Y = Y, previousObject = previousObject, newObject = newObject, 
		m.nodeMinSize = nodesize, m.maxNodes = maxnodes, m.depth = depth, m.depthControl = depthcontrol, 
		m.features = mtry, m.bootstrap = replace, m.overSampling = oversampling, 
		m.outputPerturbation = outputperturbationsampling, m.targetClass = targetclass, 
		m.rebalancedSampling = rebalancedsampling, m.regression = regression, m.ntree = ntree, m.use.OOB = OOB, m.BreimanBounds = BreimanBounds, m.treeSubsampleRate = subsamplerate, m.variableImportance = importance, 
		m.treeBagging = bagging, m.randomCombination = randomcombination, m.randomFeature = randomfeature, m.whichCatVariables = categoricalvariablesidx, m.featureSelectionRule = featureselectionrule[1], m.logX = logX, m.classwt = classwt, m.classcutoff = classcutoff, m.threads = threads, m.useSubTrees = usesubtrees,m.parallelPackage = parallelpackage[1])
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
		{ 	rebalancedsamplingString = paste(rebalancedsamplingString, rebalancedsampling[i], sep=",") }	
	}
	if (RUF.model$regression)
	{	classcutoff = c(0,0)  }	
	if (as.numeric(classcutoff[2]) == 0)
	{	classcutoffString = FALSE }
	else
	{
		classcutoffString = levels(Y)[which(levels(Y) == as.character(classcutoff[1]))]
		classcutoffString = paste("Class ", classcutoffString, "," , as.numeric(classcutoff[2])*100, "%", sep ="")
	}	
	if (is.null(previousObject))
	{
		previousObject = filter.object(previousObject)
		newObject = filter.object(newObject)		
		if (rownames(previousObject$forestParams)[1] == "reduceDimension" )
		{		
			params.big = as.character(c(as.logical(c(previousObject$forestParams[[1]][1], previousObject$forestParams[[1]][2])),previousObject$forestParams[[1]][2]))
			names.big = c("reduceDimension", "reduceAll", "nforest")
			value.big = 3 
		}
		else {  params.big = names.big = NULL; value.big = 0 }		
		paramsObject = c(params.big, ntree, mtry, nodesize, maxnodes, as.character(replace), bagging, depth, 
		ifelse(is.null(depthcontrol), length(depthcontrol)> 0, depthcontrol), OOB, importance, subsamplerate, 
		ifelse(is.null(classwt), length(classwt) > 0, classwtString), classcutoffString,
		ifelse((oversampling == 0),(oversampling != 0), oversampling), outputperturbationsampling, 
		ifelse(length(targetclass) > 1, targetclassString, targetclass), 
		ifelse(length(rebalancedsampling) > 1, rebalancedsamplingString, rebalancedsampling), 
		ifelse((randomcombination[1] == 0),(randomcombination != 0), randomcombinationString), randomfeature,  
		ifelse(is.null(categoricalvariablesidx), length(categoricalvariablesidx) > 0, length(categoricalvariablesidx)), featureselectionrule[1])								
		if (RUF.model$regression) 
		{ 
			paramsObject[21 + value.big] = "Sum of squared residuals"
			if (paramsObject[5 + value.big] == "FALSE") { paramsObject[9 + value.big] = FALSE }
		}
		names(paramsObject) = c(names.big, "ntree", "mtry", "nodesize", "maxnodes", "replace", "bagging", "depth", "depthcontrol", "OOB", "importance", "subsamplerate", "classwt", "classcutoff", "oversampling", "outputperturbationsampling", "targetclass", "rebalancedsampling", "randomcombination", "randomfeature", "categorical variables", "featureselectionrule")	
		paramsObject = as.data.frame(paramsObject)			
		RUFObject <- genericOutput(xtest, ytest, paramsObject, RUF.model, ytrain = Y, classcutoff = classcutoff)		
		if (is.null(Y) &  !RUF.model$regression)
		{  RUFObject$classes =  previousObject$classes  }		
		RUF.model$logX = logX
		if (!is.null(newObject$y)) { RUFObject$y = newObject$y	}			
		RUFObject$variablesNames = colnames(X)		
		if (!is.null(categoricalvariablesidx ))
		{  RUFObject$categoricalvariables = categoricalvariablesidx } 	
	}
	else
	{ 	
		paramsObject = cbind(previousObject$forestParams, newObject$forestParams) 
		RUF.model$logX = RUF.model$logX[1]		
		RUFObject <- genericOutput(xtest, ytest, paramsObject, RUF.model, ytrain = Y, classcutoff = classcutoff)		
		if (is.null(Y) & !RUF.model$regression) {  RUFObject$classes = previousObject$classes  }		
		RUFObject$logX = RUF.model$logX
		RUFObject$y = c(previousObject$y,newObject$y)		
		RUFObject$variablesNames = previousObject$variablesNames
	}	
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
	RUFObject$call <- match.call()
	class(RUFObject) <- "randomUniformForest"	
	RUFObject
}
		
randomUniformForestCore.merge <- function(X = NULL,Y = NULL, previousObject = NULL, newObject = NULL, 
	m.nodeMinSize = 1,
	m.maxNodes = Inf,
	m.depth = Inf,
    m.depthControl = NULL,	
	m.features = floor(4/3*ncol(X)), 
	m.bootstrap = TRUE,
	m.overSampling = 0,
	m.outputPerturbation = FALSE,
	m.targetClass = 0,
	m.rebalancedSampling = FALSE,
	m.regression = TRUE, 
	m.ntree = 100, 
	m.use.OOB = TRUE,
	m.BreimanBounds = ifelse(m.use.OOB, TRUE, FALSE),
	m.treeSubsampleRate = 1,
	m.variableImportance = TRUE,
	m.treeBagging = FALSE,
	m.classwt = NULL,	
	m.randomCombination = 0,
	m.randomFeature = FALSE,
	m.whichCatVariables = NULL,
	m.featureSelectionRule = c("entropy", "gini", "random", "L1", "L2"),
	m.logX = FALSE,
	m.classcutoff = c(0,0),
	m.useSubTrees = FALSE, 
	m.threads = "auto",
	m.parallelPackage = "doParallel")
{
    if (is.null(previousObject)) { stop ("no object to merge") }	
	if (is.null(newObject))
	{
		previousParams = previousObject$forestParams		
		if (previousParams["importance" , 1] !=  m.variableImportance)	
		{  m.variableImportance = previousParams["importance", 1] }
		if (previousParams["OOB", 1] !=  m.use.OOB)	{   m.use.OOB  = previousParams["OOB", 1]  }		
		newObject <- randomUniformForestCore(X, trainLabels = Y, depth = m.depth, depthControl = m.depthControl, 
		ntree = m.ntree, nodeMinSize = m.nodeMinSize, maxNodes = m.maxNodes, features = m.features, 
		rf.bootstrap = m.bootstrap, use.OOB = m.use.OOB, BreimanBounds = m.BreimanBounds, threads = m.threads, rf.treeSubsampleRate = m.treeSubsampleRate, rf.overSampling = m.overSampling, rf.targetClass = m.targetClass, rf.rebalancedSampling = m.rebalancedSampling, variableImportance = m.variableImportance, 
		rf.treeBagging = m.treeBagging, rf.regression = m.regression, rf.randomCombination = m.randomCombination, 
		rf.randomFeature = m.randomFeature, whichCatVariables = m.whichCatVariables, 			rf.outputPerturbationSampling = m.outputPerturbation, logX = m.logX, classwt = m.classwt, 
		classCutOff = m.classcutoff, rf.featureSelectionRule = m.featureSelectionRule, useSubTrees = m.useSubTrees,
		parallelPackage = m.parallelPackage[1])
	}			
	previousObject = filter.forest(previousObject)
	newObject = filter.forest(newObject)	
	if (length(previousObject) != length(newObject))
	{	
		objectNames = c("OOB.strengthCorr", "OOB.votes", "OOB.predicts", "OOB",  "pred.error", "percent.varExplained")
		if (is.null(previousObject$variableImportance) | is.null(newObject$variableImportance))
		{	objectNames = c(objectNames, "variableImportance") }		
		cat("'OOB object' and\\or 'Breiman Bounds' of randomUniformForest are not present in all models,\n hence they both will be dropped in the combined model. To keep them, please enable 'OOB' (and disable 'Breiman Bounds') option\n in all your models, or enable (or disable) both options in all your models.\n")		
		previousObject = rmInAListByNames(previousObject, objectNames)
		newObject = rmInAListByNames(newObject,objectNames)
	}	
	return( onlineCombineRUF(previousObject, newObject) )
}

rUniformForest.combine <- function(...)
{
	object <- list(...)
	if  (length(object) == 1)	{	object = object[[1]]	}		
	n <- length(object)	
	if (n > 1)	{	bigObject <- rUniformForest.merge(previousObject = object[[1]], newObject = object[[2]])	}
	else 	{ stop ("combine need at least two random uniform forests objects") }	
	if (n > 2)
	{
		for (i in 3:n) { bigObject <- rUniformForest.merge(previousObject = bigObject, newObject = object[[i]]) }
	}
	return(bigObject)
}


rUniformForest.grow <- function(object, X, Y = NULL,  ntree = 100, threads = "auto")
{
	if (ntree < 10)	{ stop ("please provide at least 10 trees to add") }	
	object <- filter.object(object)
	paramsObject <- object$forestParams	
	if (is.null(Y))	{  Y = object$y[1:nrow(X)]  }	
	if (!is.null(object$predictionObject))
	{	
		object = rmInAListByNames(object, c("errorObject", "predictionObject"))
		class(object) = "randomUniformForest"	
	}	
	callArgs <- as.character(object$call)	
	if ((!object$forest$regression) & !is.factor(Y)) { Y = as.factor(Y) }
	BreimanBounds = TRUE
	if (is.null(object$forest$OOB.strengthCorr)) {	BreimanBounds = FALSE	}		
	if (length(callArgs) <= 3)
	{ 
		if ( !is.null(object$formula))
		{	object2 <- randomUniformForest.formula(object$formula, data = X, ntree = ntree, BreimanBounds = BreimanBounds)	}
		else
		{	object2 <- randomUniformForest.default(X, Y = Y, ntree = ntree, BreimanBounds = BreimanBounds) 	}
	}
	else
	{
		usesubtrees = if (length(which(callArgs == "usesubtrees")) > 0) { TRUE } else { FALSE }
		bagging = if (length(which(callArgs == "bagging")) > 0) { TRUE } else { FALSE }
		na.action = { if (length(which(callArgs == "omit")) > 0) { "omit" } 
					else 
					{ 
						if (length(which(callArgs == "accurateImpute")))
						{ "accurateImpute" }
						else
						{ "fastImpute"}
					} }
		
		classCutOffString  = as.character(paramsObject["classcutoff", 1])
		classcutoff = rep(0,2)
	    if (classCutOffString != "FALSE")
		{
			classCutOffString = rm.string(classCutOffString, "%")
			classCutOffString = rm.string(classCutOffString, "Class ")
			classCutOffString = strsplit(classCutOffString,",")
			classcutoff[1] = which(object$classes == classCutOffString[[1]][1])
			classcutoff[2] = as.numeric(classCutOffString[[1]][2])/100
		}		
		if ( !is.null(object$formula))
		{
			object2 <- randomUniformForest.formula(object$formula, data = X, ntree = ntree, 
				mtry = ifelse(bagging, ncol(X), as.numeric(as.character(paramsObject["mtry",1]))),
				usesubtrees = ifelse(usesubtrees, TRUE, FALSE),
				nodesize = as.numeric(as.character(paramsObject["nodesize",1])), 
				maxnodes = as.numeric(as.character(paramsObject["maxnodes",1])),
				depth = as.numeric(as.character(paramsObject["depth",1])),
				depthcontrol = if (is.logical(as.logical(as.character(paramsObject["depthcontrol",1])))) {  NULL } 
					else 
					{ 
						if (is.numeric(as.character(paramsObject["depthcontrol",1]))) 
						{	 as.numeric(as.character(paramsObject["depthcontrol",1])) }
						else { "random" }
					},
				regression = object$forest$regression,
				replace = as.logical(as.character(paramsObject["replace",1])),
				OOB = as.logical(as.character(paramsObject["OOB",1])),
				BreimanBounds = BreimanBounds,
				subsamplerate = as.numeric(as.character(paramsObject["subsamplerate",1])),
				importance = as.logical(as.character(paramsObject["importance",1])),
				bagging = ifelse(bagging, TRUE, FALSE),
				unsupervised = if (is.null(object$unsupervised)) { FALSE } else { TRUE }, 
				unsupervisedMethod = if (is.null(object$unsupervised)) { unsupervisedMethod = c("uniform univariate sampling", "uniform multivariate sampling", "with bootstrap") } else { object$unsupervisedMethod },
				classwt = if (as.character(paramsObject["classwt", 1]) == "FALSE") { NULL } 
					else {  as.numeric(strsplit(as.character(paramsObject["classwt", 1]),",")[[1]]) },
				oversampling = if (as.character(paramsObject["oversampling", 1]) == "FALSE") { 0 } 
					else { as.numeric(as.character(paramsObject["oversampling", 1])) },
				targetclass = as.numeric(as.character(paramsObject["targetclass", 1])),
				outputperturbationsampling = if (as.character(paramsObject["outputperturbationsampling", 1]) == "FALSE") { FALSE } 	else { TRUE },
				rebalancedsampling = if (as.character(paramsObject["rebalancedsampling", 1]) == "FALSE") 
					{ as.logical(as.character(paramsObject["rebalancedsampling", 1])) } 
					else 
					{ 
						if (length(paramsObject["rebalancedsampling", 1]) > 1)	{  as.numeric(as.character(paramsObject["rebalancedsampling", 1])) }
						else {  TRUE }
					},
				featureselectionrule = if (object$forest$regression) 
									   {	
											if (as.character(paramsObject["featureselectionrule", 1]) != "Sum of squared residuals" ) 
											{  "L1"  }
											else
											{ "L2"  }
									   }
									   else
									   {	as.character(paramsObject["featureselectionrule", 1])	},			
				randomcombination = if (as.character(paramsObject["randomcombination", 1]) == "FALSE") { 0 } 
					else {	as.numeric(strsplit(as.character(paramsObject["randomcombination", 1]),",")[[1]]) },
				randomfeature = if (as.character(paramsObject["randomfeature", 1]) == "FALSE") { FALSE } 
					else { TRUE },
				categoricalvariablesidx = if (as.character(paramsObject["categorical variables", 1]) == "FALSE") { NULL } 
					else { object$categoricalvariables },
				na.action = na.action,
				logX = object$forest$logX,			
				classcutoff = classcutoff,
				threads = threads,
				parallelpackage = "doParallel"
			)		
		}
		else
		{
			object2 <- randomUniformForest.default(X, Y = Y, ntree = ntree, 
				mtry = ifelse(bagging, ncol(X), as.numeric(as.character(paramsObject["mtry",1]))),
				usesubtrees = ifelse(usesubtrees, TRUE, FALSE),
				nodesize = as.numeric(as.character(paramsObject["nodesize",1])),
				maxnodes = as.numeric(as.character(paramsObject["maxnodes",1])),
				depth = as.numeric(as.character(paramsObject["depth",1])),
				depthcontrol = 	if (is.logical(as.logical(as.character(paramsObject["depthcontrol",1])))) 
								{  NULL } 
								else 
								{ 
									if (is.numeric(as.character(paramsObject["depthcontrol",1])))
									{  as.numeric(as.character(paramsObject["depthcontrol",1])) }
									else 
									{ "random" }
								},
				regression = object$forest$regression,
				replace = as.logical(as.character(paramsObject["replace",1])),
				OOB = as.logical(as.character(paramsObject["OOB",1])),
				BreimanBounds = BreimanBounds,
				subsamplerate = as.numeric(as.character(paramsObject["subsamplerate",1])),
				importance = as.logical(as.character(paramsObject["importance",1])),
				bagging = ifelse(bagging, TRUE, FALSE),
				unsupervised = if (is.null(object$unsupervised)) { FALSE } else { TRUE }, 
				unsupervisedMethod = if (is.null(object$unsupervised)) { unsupervisedMethod = c("uniform univariate sampling", "uniform multivariate sampling", "with bootstrap") } else { object$unsupervisedMethod },
				classwt = if (as.character(paramsObject["classwt", 1]) == "FALSE") { NULL } 
						else {  as.numeric(strsplit(as.character(paramsObject["classwt", 1]),",")[[1]]) },
				oversampling = if (as.character(paramsObject["oversampling", 1]) == "FALSE") { 0 } 
							else { as.numeric(as.character(paramsObject["oversampling", 1])) },
				targetclass = as.numeric(as.character(paramsObject["targetclass", 1])),
				outputperturbationsampling = if (as.character(paramsObject["outputperturbationsampling", 1]) == "FALSE") { FALSE } 
											 else { TRUE },
				rebalancedsampling = if (as.character(paramsObject["rebalancedsampling", 1]) == "FALSE") 
									{ as.logical(as.character(paramsObject["rebalancedsampling", 1])) } 
									else 
									{ 
										if (length(paramsObject["rebalancedsampling", 1]) > 1)	
										{  as.numeric(as.character(paramsObject["rebalancedsampling", 1])) }
										else 
										{  TRUE }
									},
				featureselectionrule = 	if (object$forest$regression) 
										{	
											if (as.character(paramsObject["featureselectionrule", 1]) != "Sum of squared residuals" ) 
											{  "L1"  }
											else
											{ "L2"  }
										}
										else
										{	as.character(paramsObject["featureselectionrule", 1])	},		
				randomcombination = if (as.character(paramsObject["randomcombination", 1]) == "FALSE") { 0 } 
									else {	as.numeric(strsplit(as.character(paramsObject["randomcombination", 1]),",")[[1]]) },
				randomfeature = if (as.character(paramsObject["randomfeature", 1]) == "FALSE") { FALSE } 
								else { TRUE },
				categoricalvariablesidx = if (as.character(paramsObject["categorical variables", 1]) == "FALSE") { NULL } 
										  else { object$categoricalvariables },
				na.action = na.action,
				logX = object$forest$logX,			
				classcutoff = classcutoff,
				threads = threads,
				parallelpackage = "doParallel"
			)
		}
	}	
	objectCombined <- rUniformForest.combine(object, object2)
	objectCombined$call = object$call
	cat("OOB evaluation might become obsolete if training sample remains unchanged.\n")
	return(objectCombined)
}

rm.trees <- function(rufObject, X = NULL, Y = NULL, 
	method = c("default", "random", "oldest", "newest", "optimal", "quantile"), 
	howMany = NULL, 
	rm.sample = 0.1)
{	
	tmp.rufObject = NULL
	if (!is.null(rufObject$time.elapse)) 
	{	
		tmp.rufObject = rufObject
		rufObject = filter.object(rufObject)	
	}
	ntree = if (!is.null(rufObject$object)) { length(rufObject$object$forest$object) } else { length(rufObject$forest$object) }
	regression = rufObject$forest$regression
	if (method[1] == "default")
	{			
		if (is.null(rufObject$forest$OOB.votes))
		{ stop("Forest can not be optimized without OOB data") }		
		treesAvgScore <- apply(rufObject$forest$OOB.votes, 2, function(Z) length(rmNA(match(Z, rufObject$y))))
		cuttingEdge <- round(quantile(treesAvgScore, rm.sample),0)		
		idx = which(treesAvgScore <= cuttingEdge)
		idx2 = c(idx,2*idx)		
		rufObject$forest$object = rm.InAList(rufObject$forest$object, idx2)
	}		
	if (method[1] == "random")
	{			
		rm.idx = floor(rm.sample*ntree)
		idx = sample(ntree,rm.idx)
		
		rufObject$forest$object = rm.InAList(rufObject$forest$object,idx)
	}	
	if (method[1] == "oldest")
	{			
		if (is.null(howMany)) { stop("No number of trees has been found") }
		if (!is.numeric(howMany)) { stop("number of trees must be numeric") }
		rufObject$forest$object = rm.InAList(rufObject$forest$object, (ntree - howMany + 1):ntree)
	}	
	if (method[1] == "newest")
	{			
		if (is.null(howMany)) { stop("No number of trees has been found") }
		if (!is.numeric(howMany)) { stop("number of trees must be numeric") }
		rufObject$forest$object = rm.InAList(rufObject$forest$object, 1:(ntree - howMany))
	}	
	if (method[1] == "quantile")
	{	
		ntreeLength = unlist(lapply(rufObject$forest$object, nrow))
		idx = which(ntreeLength >= quantile(ntreeLength, 1 - rm.sample/2) |  ntreeLength <= quantile(ntreeLength, rm.sample/2) )		
		rufObject$forest$object = rm.InAList(rufObject$forest$object,idx)
	}	
	if (method[1] == "optimal")
	{	
		if (!is.null(rufObject$OOB.votes)) 	{  OOBVotes = rufObject$OOB.votes	}
		else 
		{	
			if (is.null(X)) { stop( "there is no data for optimal forest") }
			else {  predictRF = predict(rufObject, X, type = "all"); OOBVotes = predictRF$all.votes;  OOBPredicts = predictRF$majority.vote  }
		}
		if (regression)
		{ 
			L2error = mean(rufObject$forest$OOB.strengthCorr$PE.tree, na.rm =TRUE)						
			if (!is.null(rufObject$forest$object$OOB.votes)) { OOBVotes[OOBVotes == Inf] = NA }
			L2errorTree = colMeans(abs(OOBVotes), na.rm = TRUE)
			idx = which(L2errorTree > L2error)
			if (length(idx) == ntree) {  idx = which(L2errorTree  > mean(L2errorTree) ) }
		}
		else
		{ 
			TreeImportance = localTreeImportance(rufObject, OOBVotes = OOBVotes, OOBPredicts = OOBPredicts)$pweights
		    idx = which(TreeImportance < mean(TreeImportance))
		}
		rufObject$forest$object = rm.InAList(rufObject$forest$object, idx)
	}
	cat("Trees have been removed, but old number of trees will be printed\n")
	if (is.null(tmp.rufObject))	{	return(rufObject) 	}
	else
	{   
		tmp.rufObject$object = rufObject
		return(tmp.rufObject)  
	}
}
#END OF FILE