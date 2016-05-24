.womComputeDIC=function(likeVec,likePostMean){
  Dbar = -2*mean(likeVec)
  Dthetabar = -2*likePostMean
  pD = Dbar - Dthetabar
  DIC = pD + Dbar
  x=array(c(DIC,pD))
  names(x)=c("DIC","pD")
  return(x)
}


womRlogPosterior = function(x,setup) 
{
	.Call("RLogPosteriorNoCov",x,
		as.integer(setup$model$newDat2Cat[,1]),
		as.integer(setup$model$newDat2Cat[,2]),
		as.integer(setup$model$newDat2Cat[,3]),
		as.integer(setup$model$newDat2Cat[,4]),
		as.integer(as.character(setup$model$newDat2Cat[,5])),
		.womIntegerMatrix(setup$model$newDat2Cat[,-(1:5)]),
		as.matrix(setup$model$newDat2Cont[,-(1:5)]),
		as.integer(setup$model$effects[1,]),
		as.integer(setup$model$incCont),
		as.integer(setup$model$effects[1,]*0),
		as.integer(setup$model$effects[2,]),
		as.integer(setup$model$incCont),
		as.integer(setup$model$effects[2,]*0),
		as.integer(setup$model$effects[3,]),
		as.integer(setup$model$incCont),
		as.integer(setup$model$effects[3,]*0),
		setup$priors$IGa0,
		setup$priors$IGb0,
		setup$priors$muKMean,
		setup$priors$muKSD^2,
		setup$priors$muZMean,
		setup$priors$muZSD^2,
		setup$priors$muGMean,
		setup$priors$muGSD^2, 1,
		as.integer(setup$Ktype),package="WMCapacity") 
}

womRgradLogPosterior = function(x,setup)
{
	.Call("RgradLogPosteriorNoCov",x,
		as.integer(setup$model$newDat2Cat[,1]),
		as.integer(setup$model$newDat2Cat[,2]),
		as.integer(setup$model$newDat2Cat[,3]),
		as.integer(setup$model$newDat2Cat[,4]),
		as.integer(as.character(setup$model$newDat2Cat[,5])),
		.womIntegerMatrix(setup$model$newDat2Cat[,-(1:5)]),
		as.matrix(setup$model$newDat2Cont[,-(1:5)]),
		as.integer(setup$model$effects[1,]),
		as.integer(setup$model$incCont),
		as.integer(setup$model$effects[1,]*0),
		as.integer(setup$model$effects[2,]),
		as.integer(setup$model$incCont),
		as.integer(setup$model$effects[2,]*0),
		as.integer(setup$model$effects[3,]),
		as.integer(setup$model$incCont),
		as.integer(setup$model$effects[3,]*0),
		setup$priors$IGa0,
		setup$priors$IGb0,
		setup$priors$muKMean,
		setup$priors$muKSD^2,
		setup$priors$muZMean,
		setup$priors$muZSD^2,
		setup$priors$muGMean,
		setup$priors$muGSD^2, 1,
		as.integer(setup$Ktype),package="WMCapacity") 
}


womRlogLikelihood = function(x,setup){ 
.Call("RLogLikelihood",x, as.integer(setup$model$newDat2Cat[,1]),
      as.integer(setup$model$newDat2Cat[,2]),
      as.integer(setup$model$newDat2Cat[,3]),
      as.integer(setup$model$newDat2Cat[,4]),
      as.integer(as.character(setup$model$newDat2Cat[,5])),
      .womIntegerMatrix(setup$model$newDat2Cat[,-(1:5)]),
      as.matrix(setup$model$newDat2Cont[,-(1:5)]),
      as.integer(setup$model$effects[1,]),
      as.integer(setup$model$effects[2,]),
      as.integer(setup$model$effects[3,]), 1,
      as.integer(setup$Ktype),package="WMCapacity") }

womRPredVals = function(x,setup){ 
.Call("RPredictedProbabilities", x,
      as.integer(setup$model$newDat2Cat[,1]),
      as.integer(setup$model$newDat2Cat[,2]),
      as.integer(setup$model$newDat2Cat[,3]),
      as.integer(setup$model$newDat2Cat[,4]),
      as.integer(as.character(setup$model$newDat2Cat[,5])),
      .womIntegerMatrix(setup$model$newDat2Cat[,-(1:5)]),
      as.matrix(setup$model$newDat2Cont[,-(1:5)]),
      as.integer(setup$model$effects[1,]),
      as.integer(setup$model$effects[2,]),
      as.integer(setup$model$effects[3,]), 1,
      as.integer(setup$Ktype),package="WMCapacity") }


womRlogPosteriorWithCov = function(x,setup,precList,means) 
{
	.Call("RLogPosteriorWithCov",x,
		as.integer(setup$model$newDat2Cat[,1]),
		as.integer(setup$model$newDat2Cat[,2]),
		as.integer(setup$model$newDat2Cat[,3]),
		as.integer(setup$model$newDat2Cat[,4]),
		as.integer(as.character(setup$model$newDat2Cat[,5])),
		.womIntegerMatrix(setup$model$newDat2Cat[,-(1:5)]),
		as.matrix(setup$model$newDat2Cont[,-(1:5)]),
		as.integer(setup$model$effects[1,]),
		as.integer(setup$model$incCont),
		as.integer(setup$model$effects[1,]*0),
		as.integer(setup$model$effects[2,]),
		as.integer(setup$model$incCont),
		as.integer(setup$model$effects[2,]*0),
		as.integer(setup$model$effects[3,]),
		as.integer(setup$model$incCont),
		as.integer(setup$model$effects[3,]*0),
		setup$priors$IGa0,
		setup$priors$IGb0,
		setup$priors$muKMean,
		setup$priors$muKSD^2,
		setup$priors$muZMean,
		setup$priors$muZSD^2,
		setup$priors$muGMean,
		setup$priors$muGSD^2, 1,
		as.integer(setup$Ktype),
		precList, means,
		as.integer(setup$model$covObs),
		as.integer(setup$model$covSizes),
		as.integer(setup$model$parStart),
		as.integer(setup$model$effSlope), package="WMCapacity") 
}

womRgradLogPosteriorWithCov = function(x,setup,precList,means)
{
	.Call("RgradLogPosteriorWithCov",x,
		as.integer(setup$model$newDat2Cat[,1]),
		as.integer(setup$model$newDat2Cat[,2]),
		as.integer(setup$model$newDat2Cat[,3]),
		as.integer(setup$model$newDat2Cat[,4]),
		as.integer(as.character(setup$model$newDat2Cat[,5])),
		.womIntegerMatrix(setup$model$newDat2Cat[,-(1:5)]),
		as.matrix(setup$model$newDat2Cont[,-(1:5)]),
		as.integer(setup$model$effects[1,]),
		as.integer(setup$model$incCont),
		as.integer(setup$model$effects[1,]*0),
		as.integer(setup$model$effects[2,]),
		as.integer(setup$model$incCont),
		as.integer(setup$model$effects[2,]*0),
		as.integer(setup$model$effects[3,]),
		as.integer(setup$model$incCont),
		as.integer(setup$model$effects[3,]*0),
		setup$priors$IGa0,
		setup$priors$IGb0,
		setup$priors$muKMean,
		setup$priors$muKSD^2,
		setup$priors$muZMean,
		setup$priors$muZSD^2,
		setup$priors$muGMean,
		setup$priors$muGSD^2, 1,
		as.integer(setup$Ktype),
		precList, means,
		as.integer(setup$model$covObs),
		as.integer(setup$model$covSizes),
		as.integer(setup$model$parStart),
		as.integer(setup$model$effSlope), package="WMCapacity") 
}



