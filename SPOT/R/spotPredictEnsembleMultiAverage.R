###################################################################################
#' spotEnsembleMultiAverage
#' 
#' Weighted Average of Model Prediction \cr
#' This "multi ensemble" (see details) combines the response of the models in the ensemble using a weighted average. 
#' The predicted values of each model are weighted by the models feedback (e.g. error), and then combined by averaging.\cr
#' Relevant configuration parameters of this ensemble are: \cr
#' The function to calculate feedback (e.g. the error) of the models, default is \code{spotConfig$seq.ensemble.feed.func<-spotFeedback.error.full} \cr
#' The function average model predictions, default is \code{spotConfig$seq.ensemble.average.func<-mean} \cr
#' The number of designs created for validation, which can be 0 or larger integers, e.g \code{spotConfig$seq.ensemble.cut.num <- 0} 
#' The the cut-off type used for validation, which can be 0 or larger integers, e.g \code{spotConfig$seq.ensemble.cut.num <- "perc"} 
#' The size of the point-set cut off for validation, which can be 0 or larger integers, e.g \code{spotConfig$seq.ensemble.cut.num <- 10} 
#'
#' This is a "multi ensemble", meaning that in every sequential step all models in the ensemble are trained and evaluated.
#' The target is to actively combine all models responses, to get the best estimate on which candidate points are optimal.\cr
#' The models used are specified in the \code{spotConfig} list, for instance:\cr
#' \code{spotConfig$seq.ensemble.predictors = c(spotPredictRandomForest, spotPredictEarth, spotPredictForrester, spotPredictDace)}\cr
#' To specify the settings of each individual model, set:\cr
#' \code{seq.ensemble.settings = list(list(setting=1),list(setting=2),list(setting=3),list(setting=4))}\cr
#' Any parameters set in each of the corresponding lists (here: 4 individual lists) will overwrite settings in the main \code{spotConfig} list,
#' when the concerned model function is called.
#'
#' @param rawB unmerged data
#' @param mergedB merged data
#' @param design new design points which should be predicted
#' @param spotConfig global list of all options, needed to provide data for calling functions
#' @param fit if an existing model ensemble fit is supplied, the models will not be build based on 
#'				data, but only evaluated with the existing fits (on the design data). To build the model, 
#'				this parameter has to be NULL. If it is not NULL the parameters mergedB and rawB will not be 
#'				used at all in the function.
#'
#' @return returns the list \code{spotConfig} 
#'
#' @references - M. Friese, M. Zaefferer, T. Bartz-Beielstein, O. Flasch, P. Koch, W. Konen, and
#' B. Naujoks. Ensemble based optimization and tuning algorithms. In F. Hoffmann
#' and E. Huellermeier, editors, Proceedings 21. Workshop Computational Intelligence,
#' p. 119-134. Universitaetsverlag Karlsruhe, 2011.
#'
#' @export
###################################################################################
spotEnsembleMultiAverage <- function(rawB,mergedB,design,spotConfig,fit=NULL){	
	design <- spotInitializePredictor(design,"data.frame",spotConfig$alg.roi,NULL,"spotEnsembleMultiAverage",spotConfig$io.verbosity)
	# Initialize
	if(is.null(spotConfig$seq.ensemble.feed.func)){spotConfig$seq.ensemble.feed.func<-spotFeedback.error.full}
	if(is.null(spotConfig$seq.ensemble.average.func)){spotConfig$seq.ensemble.average.func<-mean}
	
	if(is.null(fit)){ #compute errors and build models if no existing information is provided

		# Compute Error from last steps models
		if(!is.null(spotConfig$seq.modelFit$fitList)){ #can only be done if fitlist exists from last step
			err<-spotConfig$seq.ensemble.feed.func(spotConfig,mergedB,rawB)
		}
		else{
			err=NA;
		}
		#Build Models:
		res<-spotEnsembleModelBuilding(rawB, mergedB, design, spotConfig) 

	}else{ #Else do only prediction on preexisting information
		err <- spotConfig$seq.modelFit$err
		res <- spotEnsembleModelEvaluate(spotConfig,design,all=TRUE)	
	}

	#Check error
	if(is.na(err[1])){ #if no error from last step, calculate from this step (i.e. first step)
		spotConfig$seq.modelFit$fitList<-res$fits
		err=spotConfig$seq.ensemble.feed.func(spotConfig,mergedB,rawB)
	}
	
	results <- as.data.frame(res$Y)
	
	#Calculate weighted output of all models
	nn=ncol(results)
	nmodels=length(spotConfig$seq.ensemble.predictors)
	if(nn>nmodels){weights=rep(err,each=nn/nmodels)} #replicate weights
	else{weights=err} #no replicate needed	
	results <- t(t(results)*weights) #weighted results

	spotWriteLines(spotConfig$io.verbosity,3,"spotEnsembleMultiAverage finished");

	# Determine values to be returned inside spotConfig
	spotConfig$seq.largeDesignY<-as.data.frame(apply(results,1,spotConfig$seq.ensemble.average.func));
	spotConfig$seq.modelFit$err<-err
	spotConfig$seq.modelFit$fitList<-res$fits
	spotConfig
}

