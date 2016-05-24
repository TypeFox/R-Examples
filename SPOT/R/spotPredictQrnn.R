###################################################################################
#' Meta Model Interface: Quantile Regression Neural Network
#' 
#' This meta model uses the "qrnn" package to build a quantile regression neural network.
#'
#' @param rawB unmerged data
#' @param mergedB merged data
#' @param design new design points which should be predicted
#' @param spotConfig global list of all options, needed to provide data for calling functions
#' @param fit if an existing model fit is supplied, the model will not be build based on 
#'				data, but only evaluated with the model fit (on the design data). To build the model, 
#'				this parameter has to be NULL. If it is not NULL the parameters mergedB and rawB will not be 
#'				used at all in the function.
#'
#' @return returns the list \code{spotConfig} with two new entries:\cr
#' 	spotConfig$seq.modelFit fit of the earth model used with predict() \cr
#'	spotConfig$seq.largeDesignY the y values of the design, evaluated with the fit
#' @export
###################################################################################
spotPredictQrnn<- function(rawB,mergedB,design,spotConfig,fit=NULL){	
	design <- spotInitializePredictor(design,"as.matrix",spotConfig$alg.roi,"qrnn","spotPredictQrnn",spotConfig$io.verbosity)	
	########################################################
	# BUILD
	########################################################	
	if(is.null(fit)){
		xNames <- row.names(spotConfig$alg.roi); 	
		yNames <- spotConfig$alg.resultColumn
		x <- rawB[xNames]		
		if(length(yNames)==1){
			y <- rawB[[yNames]]
			fit <- qrnn::qrnn.fit(as.matrix(x),as.matrix(y),n.hidden=5)		
		}
		else{#Distinction for multi criteria spot 
			y <- rawB[yNames]
			fit=list()
			res=list()
			for (i in 1:length(yNames)){
				fit[[i]]<-qrnn::qrnn.fit(as.matrix(x),as.matrix(y[,i]),n.hidden=5)	
			}			
		}	
	}else{
		fit<-fit
	}
	########################################################
	# PREDICT
	########################################################
	if(!is.null(design)){ 
		if(length(spotConfig$alg.resultColumn)>1){
			res=list()
			for (i in 1:length(fit)){
				res[[i]]<-qrnn::qrnn.predict(design,fit[[i]])
			}	
		}
		else{res <- qrnn::qrnn.predict(design,fit)}		
		res[which(is.na(res)==TRUE)] = median(fit$y,na.rm = TRUE) #replaces NA values with median of y..?
	}else{res <- NULL}
	########################################################
	# OUTPUT
	########################################################
	spotWriteLines(spotConfig$io.verbosity,3,"spotPredictQrnn finished")
	spotConfig$seq.modelFit<-fit
	spotConfig$seq.largeDesignY<-as.data.frame(res)
	spotConfig
}
