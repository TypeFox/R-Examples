###################################################################################
#' Meta Model Interface: Neural network
#' 
#' Prediction based on neuralnet package.
#' Can be used both for single and multi objective SPOT. 
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
#' @seealso \code{\link{spotPredictLmFactor}} 
#' @export
###################################################################################
spotPredictNeuralnet <- function(rawB,mergedB,design,spotConfig,fit=NULL){		
	design <- spotInitializePredictor(design,"data.frame",spotConfig$alg.roi,"neuralnet","spotPredictNeuralnet",spotConfig$io.verbosity)
	xNames <- row.names(spotConfig$alg.roi)
	########################################################
	# BUILD
	########################################################
	if(is.null(fit)){
		yNames <- spotConfig$alg.resultColumn
		frml = paste(paste(yNames,collapse="+"),paste(xNames,collapse="+"),sep="~")
		fit <- neuralnet::neuralnet(frml,rawB,hidden=2,threshold= 0.001,stepmax=1e5)	
		while(is.null(fit$weights)){
			fit <- neuralnet::neuralnet(frml,rawB,hidden=2,threshold= 0.001,stepmax=1e5)	
		}
	}
	########################################################
	# PREDICT
	########################################################
	if(!is.null(design)){ 	
		res <-neuralnet::compute(fit,design)$net.result
	}else{res <- NULL}	
	########################################################
	# OUTPUT
	########################################################	
	spotWriteLines(spotConfig$io.verbosity,3,"spotPredictNeuralnet finished")
	spotConfig$seq.modelFit<-fit
	spotConfig$seq.largeDesignY<-as.data.frame(res)
	spotConfig
}

