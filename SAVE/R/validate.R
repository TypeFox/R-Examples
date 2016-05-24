##########################################################################
## Validation Method
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 3, April 2013.
##
## Copyright (C) 2013-present Jesus Palomo, Gonzalo Garcia-Donato, 
##							  and Rui Paulo
##    
##########################################################################

validate.SAVE <- function(object, newdesign=NULL, calibration.value="mean", prob=0.90, n.burnin=0, n.thin=1, tol=1E-10){

	if (object@constant.controllables && (!is.null(newdesign))){
		stop("In the situation with constant controllable inputs, newdesign should be set to NULL.") 
	}
	
	if (!object@constant.controllables){
		tmp<- predictreality(object=object, newdesign=newdesign, n.burnin=n.burnin, n.thin=n.thin, tol=tol)
	}
	else{
		tmp<- predictreality(object=object, n.burnin=n.burnin, n.thin=n.thin, tol=tol)
	}
	reality<- (tmp@biaspred+tmp@modelpred)
	meanreality<- apply(X=reality, MARGIN=2, FUN=mean)
	
	# tolerance bounds
	tmpdata <- matrix(meanreality,ncol=dim(reality)[2],nrow=dim(reality)[1],
					  byrow=T)
	tmpdata2 <- reality - tmpdata
	tmpdata <- apply(X=tmpdata2,MARGIN=2,FUN=abs)
	tau.real <- apply(X=tmpdata,MARGIN=2,FUN=quantile,probs=prob)
	
	#Specify the value of the calibration parameters:
	if (length(object@calibrationnames) != 0){
		if (is.character(calibration.value)){
			if (calibration.value=="mean")
			{calibration.value<- apply(as.matrix(object@mcmcsample[,object@calibrationnames]), MARGIN=2, FUN=mean)}
			else {
				if (calibration.value=="median")
				{calibration.value<- apply(as.matrix(object@mcmcsample[,object@calibrationnames]), MARGIN=2, FUN=median)}
				else {stop("Invalid calibration.value\n")}
			}
		}
		else{
			calibration.value<- vapply(calibration.value, FUN=function(x){x[1]}, FUN.VALUE=c(0))[object@calibrationnames]
		}

		if (!object@constant.controllables){
		newdesignpure<- cbind(newdesign, matrix(calibration.value, nrow=dim(newdesign)[1], ncol=length(calibration.value), byrow=T))
		names(newdesignpure) <- c(names(newdesign),object@calibrationnames)
		}
		else{
		newdesignpure<- matrix(calibration.value, nrow=1, ncol=length(calibration.value), byrow=T)
		colnames(newdesignpure) <- object@calibrationnames
		newdesignpure<- as.data.frame(newdesignpure)
		}
	} 
	else 
		if (!object@constant.controllables){
	    newdesignpure <- newdesign
		}
		else{
		newdesignpure <- NULL
		}

	puremodel<- predictcode(object=object, newdesign=newdesignpure, n.iter=10, sampledraws=F, tol=tol)@modelmean
		
	# tolerance bounds
	tmpdata <- matrix(puremodel,ncol=dim(reality)[2],nrow=dim(reality)[1],
					  byrow=T)
	tmpdata2 <- reality - tmpdata
	tmpdata2 <- apply(X=tmpdata2,MARGIN=2,FUN=abs)
	tau.pure <- apply(X=tmpdata2,MARGIN=2,FUN=quantile,probs=prob)

	#the estimates of the bias:
	bias<- apply(X=as.matrix(reality-tmpdata), MARGIN=2, FUN=mean)
	biasL<- apply(X=as.matrix(reality-tmpdata), MARGIN=2, FUN=quantile, probs=(1-prob)/2)
	biasU<- apply(X=as.matrix(reality-tmpdata), MARGIN=2, FUN=quantile, probs=(1+prob)/2)
	
	#Results are being stored in an object called results
	result<- new("validate.SAVE")
	result@call<- object@call
	result@bayesfitcall<- object@bayesfitcall
	# To deprecate the unused parameters and to include in the 
	# call() all the default parameters not used in the call
	# to the function
	dprct <- .deprecate.parameters(call=sys.call(sys.parent(1)))
	result@validatecall <- as.call(dprct)
	if (length(object@calibrationnames) == 0){
		aux <- which(names(result@validatecall)=="calibration.value")
		#print (paste("The specified parameter calibration.value=",result@validatecall[aux]," has been removed since there are no calibration parameters in the model.",sep=''))
		result@validatecall <- result@validatecall[-aux]
	}
	
	result@validate<- cbind(meanreality, tau.real, puremodel, tau.pure, bias, biasL, biasU)
	colnames(result@validate)<- c("bias.corrected", "tau.bc", "pure.model", "tau.pm", "bias", "bias.Lower", "bias.Upper")
	if (!object@constant.controllables){
		rownames(result@validate)<- rownames(newdesign)
		result@newdesign<- newdesign		
	}
	else{rownames(result@validate)<- 1}
	
	unlink(paste0(object@wd,'/*'))
	
	return(result)
}

if(!isGeneric("validate")) {
	setGeneric(name = "validate",
			def = function(object, newdesign=NULL, calibration.value="mean",
					prob=0.90, n.burnin=0, n.thin=1, tol=1E-10,...) 
				standardGeneric("validate")
	)
}


setMethod("validate", "SAVE", 
		definition= function(object, newdesign,...) {
			validate.SAVE(object=object, newdesign=newdesign, calibration.value=calibration.value, prob=prob, n.burnin=n.burnin, n.thin=n.thin, tol=tol)
		}
)

setMethod("show","validate.SAVE",
		function(object){
		    smx <- summary.validate.SAVE(object)
			show.summary.validate.SAVE (smx) }
)
			
summary.validate.SAVE<- function(object, ...){
    result<- new ("summary.validate.SAVE")
	result@callSAVE <- object@call
	result@callbayes <- object@bayesfitcall
	result@callvalidate <- object@validatecall
	result@summaries <- object@validate
	return(result)
	# cat("\n")
	# cat("---------------\n")
	# cat("call to SAVE:\n")
	# print(object@call)
	# cat("---------------\n")
	# cat("call to bayesfit:\n")
	# print(object@bayesfitcall)
	# cat("---------------\n")
	# cat("call to validate:\n")
	# print(object@validatecall)
	# cat("---------------\n")
	# cat("Results:\n")
	# print(object@validate)
}

setMethod("summary","validate.SAVE",
	function(object){ summary.validate.SAVE (object) }
)

show.summary.validate.SAVE<- function(object){
	cat("\n")
	cat("---------------\n")
	cat("call to SAVE:\n")
	print(object@callSAVE)
	cat("---------------\n")
	cat("call to bayesfit:\n")
	print(object@callbayes)
	cat("---------------\n")
	cat("call to validate:\n")
	print(object@callvalidate)
	cat("---------------\n")
	cat("Results:\n")
	if (!is.null(object@summaries) || length(object@summaries)!=0){
	  print(object@summaries)
	  }
	}
	
setMethod("show","summary.validate.SAVE",
        function(object){show.summary.validate.SAVE (object) }
)

			
plot.validate.SAVE<- function(x, ...){
	summaries<- summary.validate.SAVE(x)@summaries
	par(mfrow=c(3,1))
	ind <- 1:dim(summaries)[1]
	
	what <- as.vector(summaries[,"pure.model"])
	delta <- as.vector(summaries[,"tau.pm"]) 
	plot(x=ind, y=what, xlab="Input points", 
		ylim=c(min(what-delta)*0.95,max(what+delta)*1.05), ylab="Pure-model", type="n")
	points(x=ind, y=what, cex=0.8)
	for (i in ind){
		lines(x=c(i,i),y=c(what[i]-delta[i],what[i]+delta[i]))
		}
		
	what <- as.vector(summaries[,"bias"])
	lo <- as.vector(summaries[,"bias.Lower"])
	up <- as.vector(summaries[,"bias.Upper"]) 
	plot(x=ind, y=what, xlab="Input points", 
		ylim=c(min(lo)*0.95,max(up)*1.05), ylab="Bias", type="n")
	points(x=ind, y=what, cex=0.8)
	for (i in ind){
		lines(x=c(i,i),y=c(lo[i],up[i]))
		}
	
	what <- as.vector(summaries[,"bias.corrected"])
	delta <- as.vector(summaries[,"tau.bc"]) 
	plot(x=ind, y=what, xlab="Input points", 
		ylim=c(min(what-delta)*0.95,max(what+delta)*1.05), ylab="Bias-corrected", type="n")
	points(x=ind, y=what, cex=0.8)
	for (i in ind){
		lines(x=c(i,i),y=c(what[i]-delta[i],what[i]+delta[i]))
		}
}

setMethod("plot",signature(x="validate.SAVE",y="missing"), 
	function(x, ...) {
	plot.validate.SAVE(x = x, ...)
}
)