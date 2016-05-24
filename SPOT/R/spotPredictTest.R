#######################################################################
#' Meta Model Interface: Test
#' 
#' This model does prediction by sampling on the real target function.
#' Supposed to be used for testing purpose only. Do not use as surrogate model.
#'
#' @param rawB unmerged data
#' @param mergedB merged data
#' @param largeDesign new design points which should be predicted
#' @param spotConfig global list of all options, needed to provide data for calling functions
#' @param fit if an existing model fit is supplied, the model will not be build based on 
#'				data, but only evaluated with the model fit (on the largeDesign data). To build the model, 
#'				this parameter has to be NULL. If it is not NULL the parameters mergedB and rawB will not be 
#'				used at all in the function.
#'
#' @return returns the list \code{spotConfig} with two new entries:\cr
#' 	spotConfig$seq.modelFit fit of the earth model used with predict() \cr
#'	spotConfig$seq.largeDesignY the y values of the large design, evaluated with the fit
#' @keywords internal
#######################################################################
spotPredictTest <- function(rawB,mergedB,largeDesign,spotConfig,fit=NULL){	
	if(spotConfig$alg.func=="spotFuncStartBranin"){fn<-spotBraninFunction}
	else if(spotConfig$alg.func=="spotFuncStartSixHump"){fn<-spotSixHumpFunction}
	else if(spotConfig$alg.func=="spotFuncStartSphere"){fn<-spotSphereFunction}
	else if(spotConfig$alg.func=="spotFuncStartRosenbrock"){fn<-spotRosenbrockFunction }
	else if(spotConfig$alg.func=="spotFuncStartRastrigin"){fn<-spotRastriginFunction}
	else if(spotConfig$alg.func=="spotFuncStartMexicanHat"){fn<-spotMexicanHatFunction}
	else {fn<-spotConfig$alg.tar.func}
	if(ncol(largeDesign)!=nrow(spotConfig$alg.roi)){ #ugly, but necessary because R is strange about converting one dimensional vectors to dataframes (confusing rows/columns)
			largeDesign<-t(data.frame(largeDesign))			
	}
	res <- apply(largeDesign,1,fn)	
	if(length(spotConfig$alg.resultColumn)>1){res<-t(res)}
	#fn2<-function(y){
	#	y + spotCalcNoise(y, spotConfig$spot.noise, spotConfig$spot.noise.type, spotConfig$spot.noise.minimum.at.value);
	#}
	#res <- apply(as.data.frame(res),1,FUN=fn2)
	spotConfig$seq.modelFit<-fn;
	spotConfig$seq.largeDesignY<-as.data.frame(res);
	spotConfig
}
