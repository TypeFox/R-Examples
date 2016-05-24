# S3 Generic methods
# model object supports: nn, nnet, rsnns, ...
nnetPredInt<-function(object, ...){
	UseMethod("nnetPredInt", object)
}

# S3 method for default
nnetPredInt.default<-function(object = NULL, xTrain, yTrain, yFit, node, wts, newData, alpha = 0.05 ,lambda = 0.5, funName = 'sigmoid', ...) {
  yPredInt = getPredInt(xTrain, yTrain, yFit, node, wts, newData, alpha = alpha, lambda = lambda, funName = funName)
  return(yPredInt)
}

# S3 method for neuralnet
nnetPredInt.nn<-function(object, xTrain, yTrain, newData, alpha = 0.05, lambda = 0.5, funName = 'sigmoid', ...) {
	# xTrain = object$covariate
	# yTrain = object$response
	yFit = object$net.result[[1]]
	wtsList = object$weights[[1]]
	m = length(wtsList)
	nodeNum = c()
	for (i in 1:m) {
		curNodeNum = dim(wtsList[[i]])[1] - 1
		nodeNum = c(nodeNum,curNodeNum)
	}
	nodeNum = c(nodeNum,1) # output layer
	wts = transWeightListToVect(wtsList,m)
	yPredInt = getPredInt(xTrain, yTrain, yFit, nodeNum, wts, newData, alpha = alpha, lambda = lambda, funName = funName)
	return(yPredInt)
}

# S3 method for nnet
nnetPredInt.nnet<-function(object, xTrain, yTrain, newData, alpha = 0.05, lambda = 0.5, funName = 'sigmoid', ...) {
	wts = object$wts
	nodeNum = object$n
	yFit = c(object$fitted.values)
	yPredInt = getPredInt(xTrain, yTrain, yFit, nodeNum, wts, newData, alpha = alpha, lambda = lambda, funName = funName)
	return(yPredInt)
}

# S3 method for RSNNS
nnetPredInt.rsnns<-function(object, xTrain, yTrain, newData, alpha = 0.05, lambda = 0.5, funName = 'sigmoid', ...) {
	# nodeNum
	nodeNum = c()
	nodeNum = c(nodeNum, object$nInputs)
	hiddenLayer = object$archParams$size
	nodeNum = c(nodeNum,hiddenLayer)
	nodeNum = c(nodeNum,object$nOutputs)
	# yFit
	yFit = object$fitted.values
	# wts : Wij + bi
	nInput = object$nInputs
	nUnit = as.numeric(extractNetInfo(object)$infoHeader$value[1])
	biasVect = extractNetInfo(object)$unitDefinitions$unitBias
	biasVect = biasVect[(nInput + 1):nUnit]
	wtsMatSNNS = weightMatrix(object)[,(nInput + 1):nUnit]
	wts = c()
	for (i in 1:(nUnit - nInput)) {
		wts = c(wts,biasVect[i])
		curNodeWts = wtsMatSNNS[,i]
		idx = which(curNodeWts != 0.0)
		wts = c(wts,curNodeWts[idx])
	}
	yPredInt = getPredInt(xTrain, yTrain, yFit, nodeNum, wts, newData, alpha = alpha, lambda = lambda, funName = funName)
	return(yPredInt)
}
