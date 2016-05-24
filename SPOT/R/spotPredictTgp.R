###################################################################################
#' Meta Model Interface: Treed Gaussian Processes
#'
#' This function implements a  model for prediction, based on Mat's tgp 
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
#' @return returns the list \code{spotConfig} with a new entry:\cr
#'  spotConfig$seq.modelFit fit of the model used with predict() \cr
#'	spotConfig$seq.largeDesignY the y values of the design, evaluated with the fit
#'
#' @seealso \code{\link{SPOT}}
#' @export
###################################################################################
spotPredictTgp <- function(rawB,mergedB,design,spotConfig,fit=NULL){
	design <- spotInitializePredictor(design,"as.matrix",spotConfig$alg.roi,"tgp","spotPredictTgp",spotConfig$io.verbosity)	
	########################################################
	# BUILD
	########################################################	
	if(is.null(fit)){ # no fit given, build model
		xNames <- row.names(spotConfig$alg.roi)
		yNames <- spotConfig$alg.resultColumn
		x <- rawB[xNames] 
		y <- rawB[yNames] 
		fit <- tgp::btgpllm(X=x,Z=y,XX=design,verb=spotConfig$io.verbosity,improv=c(1,spotConfig$seq.design.size))
		res=predict(fit,design)$ZZ.mean
		#res=fit$improv$rank;# TODO this must be tested, else use predict as above 	
	}
	########################################################
	# PREDICT
	########################################################
	if(!is.null(design)) res <- predict(fit,design)$ZZ.mean
	else res <- NULL
	########################################################
	# OUTPUT
	########################################################
	spotWriteLines(spotConfig$io.verbosity,3,"spotPredictTgp finished")	
	spotConfig$seq.modelFit=fit
	spotConfig$seq.largeDesignY<-as.data.frame(res)
	spotConfig
}