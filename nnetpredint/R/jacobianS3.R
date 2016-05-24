# S3 Generic methods
# model object supports: nn, nnet, rsnns, ...
jacobian<-function(object, ...){
	UseMethod("jacobian", object)
}

jacobian.default<-function(object = NULL, xTrain, W, B, m, nPara, funName = 'sigmoid', ...){
	jacobianMat = getJacobianMatrix(xTrain, W, B, m, nPara, funName)
	return (jacobianMat)
}

jacobian.nn<-function(object, xTrain, funName = 'sigmoid', ...){
	wtsList = object$weights[[1]]
	m = length(wtsList)
	nodeNum = c()
	for (i in 1:m) {
		curNodeNum = dim(wtsList[[i]])[1] - 1
		nodeNum = c(nodeNum,curNodeNum)
	}
	nodeNum = c(nodeNum,1) # output layer
	wts = transWeightListToVect(wtsList,m)
	nPara = length(wts)
	transWts = transWeightPara(wts, m, nodeNum)
	W = transWts$W
	B = transWts$B
	jacobianMat = getJacobianMatrix(xTrain,W,B,m,nPara,funName)
	return (jacobianMat)
}
# jacobMat = jacobian(nn,x)   # nObs * nPara

jacobian.nnet<-function(object, xTrain, funName = 'sigmoid', ...){
	wts = object$wts
	nodeNum = object$n
	m = length(nodeNum) - 1
	nPara = length(wts)
	transWts = transWeightPara(wts, m, nodeNum)
	W = transWts$W
	B = transWts$B
	jacobianMat = getJacobianMatrix(xTrain,W,B,m,nPara,funName)
	return (jacobianMat)
}
# jacobMat = jacobian(nnet,x)

jacobian.rsnns<-function(object, xTrain, funName = 'sigmoid', ...){
	#nodeNum
	nodeNum = c()
	nodeNum = c(nodeNum, object$nInputs)
	hiddenLayer = object$archParams$size
	nodeNum = c(nodeNum,hiddenLayer)
	nodeNum = c(nodeNum,object$nOutputs)
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
	m = length(nodeNum) - 1
	nPara = length(wts)
	transWts = transWeightPara(wts, m, nodeNum)
	W = transWts$W
	B = transWts$B
	jacobianMat = getJacobianMatrix(xTrain,W,B,m,nPara,funName)
	return (jacobianMat)
}

# jacobMat = jacobian(nnet,x)   # nObs * nPara
