###################################################################################
#' Meta Model Interface: Linear model with factors for SPOT
#'
#' This function uses the lm function shipped with R, while \code{\link{spotPredictLm}} uses the rsm package. 
#'
#' The linear model approach described here can incorporate parameters which are marked as FACTORS (i.e. categorical parameters) in the region of interest, see \code{\link{spotROI}}.
#' Please note that the design used to train the linear model should contain all levels of the factor variable. FACTORS are not ordered, and therefore are impossible to extrapolate on.
#' If new data is given in the \code{design} variable which contains unseen FACTOR levels, please note that this will probably create NA values in the prediction.
#' NA values might yield errors in your SPOT run, ending it prematurely. It is therefore recommended to build a initial design which contains at least one example
#' of each FACTOR level.
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
#' @seealso \code{\link{spotPredictLm}}
#' @export
###################################################################################
spotPredictLmFactor <- function(rawB,mergedB,design,spotConfig,fit=NULL){	
	design <- spotInitializePredictor(design,"data.frame",spotConfig$alg.roi,NULL,"spotPredictLmFactor",spotConfig$io.verbosity)
	xNames <- row.names(spotConfig$alg.roi)	
	########################################################
	# BUILD
	########################################################		
	if(is.null(fit)){
		yNames <- spotConfig$alg.resultColumn
		x <- rawB[xNames]		
		y <- rawB[yNames]
		tmp <- spotForceFactorX(x,spotConfig$alg.roi$type,xNames)
		x <- tmp$x #parameters are now factors if specified in roi
		spotConfig$alg.factor.levels <- tmp$levels #needed for later conversion of prediction locations	
		
		fit<-lm(y~.,data=cbind(y,x))
	}
	########################################################
	# PREDICT
	########################################################
	if(!is.null(design)){ 	
		design <- spotForceFactorDesign(design,spotConfig$alg.roi$type,xNames,spotConfig$alg.factor.levels)
		res <- predict(fit,design) #TODO: please note that factors which are _not_ known to the model will produce a prediction of NA!
	}else{res <- NULL}		
	########################################################
	# OUTPUT
	########################################################
	spotWriteLines(spotConfig$io.verbosity,3,"spotPredictLmFactor finished")
	spotConfig$seq.modelFit<-fit
	spotConfig$seq.largeDesignY<-as.data.frame(res)
	spotConfig
}
