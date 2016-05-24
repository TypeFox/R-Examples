###################################################################################
#' Meta Model Interface: Fields Kriging
#' 
#' Kriging meta model based on fields package.
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
#' 	spotConfig$seq.modelFit fit of the Krig model used with predict() \cr
#'	spotConfig$seq.largeDesignY the y values of the design, evaluated with the fit
#' @export
###################################################################################
spotPredictKrig <- function(rawB,mergedB,design,spotConfig,fit=NULL){	
	design <- spotInitializePredictor(design,"data.frame",spotConfig$alg.roi,"fields","spotPredictKrig",spotConfig$io.verbosity)	
	########################################################
	# BUILD
	########################################################	
	if(is.null(fit)){
		xNames <- row.names(spotConfig$alg.roi)
		yNames <- spotConfig$alg.resultColumn
		x <- rawB[xNames]
		y <- rawB[yNames]
		fit <- fields::Krig(x=x,Y=y)
	}else{
		fit<-fit
	}
	########################################################
	# PREDICT
	########################################################	
	if(!is.null(design)){ 	
		res <- predict(fit,design)	
		#error: predict.se(fit,lD)
		res[which(is.na(res)==TRUE)] = median(fit$y,na.rm = T) #replaces NA values with median of y..?
	}else{res <- NULL}			
	########################################################
	# OUTPUT
	########################################################		
	spotWriteLines(spotConfig$io.verbosity,3,"spotPredictKrig finished")
	spotConfig$seq.modelFit<-fit
	spotConfig$seq.largeDesignY<-as.data.frame(res)
	spotConfig
}
