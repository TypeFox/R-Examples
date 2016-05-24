###################################################################################
#' Meta Model Interface: Tree
#'  
#' A prediction model based on rpart, using a single tree model.
#' 
#' @param rawB unmerged data
#' @param mergedB merged data
#' @param design new design points which should be predicted
#' @param spotConfig global list of all options, needed to provide data for calling functions
#' @param fit if an existing model fit is supplied, the model will not be build based on 
#'				data, but only evaluated with the model fit (on the largeDesign data). To build the model, 
#'				this parameter has to be NULL. If it is not NULL the parameters mergedB and rawB will not be 
#'				used at all in the function.
#'
#' @return returns the list \code{spotConfig} with two new entries:\cr
#' 	spotConfig$seq.modelFit fit of the model used with predict() \cr
#'	spotConfig$seq.largeDesignY the y values of the design, evaluated with the fit
#'  
#' @seealso \code{\link{SPOT}}
#' @export
###################################################################################
spotPredictTree <- function(rawB,mergedB,design,spotConfig,fit=NULL){
	design <- spotInitializePredictor(design,"data.frame",spotConfig$alg.roi,"rpart","spotPredictTree",spotConfig$io.verbosity)
	########################################################
	# BUILD
	########################################################	
	if(is.null(fit)){ #build model if not given
		xNames <- row.names(spotConfig$alg.roi)
		yNames <- spotConfig$alg.resultColumn
		fit <- rpart(paste(yNames,"~",xNames), data= rawB)
	}
	########################################################
	# PREDICT
	########################################################	
	if(!is.null(design)) res <- predict(fit,design)		
	else res <- NULL
	########################################################
	# OUTPUT
	########################################################	
	spotWriteLines(spotConfig$io.verbosity,3,"spotPredictTree finished")
	spotConfig$seq.modelFit<-fit
	spotConfig$seq.largeDesignY<-as.data.frame(res)
	spotConfig
}