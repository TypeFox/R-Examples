###################################################################################
#' spotEnsembleMultiChoose
#' 
#' Choose a model based on feedback\cr
#' This "multi ensemble" (see details) evaluates each models feedback (e.g. error), 
#' and the best one is chosen to be used for the determination of new optimal candidate points. 
#' In order to allow for randomization, the epsilon parameter can be set to values between 0  and 1, where 1 means completely random (feedback is ignored)
#' and 0 means completely deterministic (choice only dependent on feedback, not random). \cr
#' Relevant configuration parameters of this ensemble are: \cr
#' The function to calculate feedback (e.g. the error) of the models, default is \code{spotConfig$seq.ensemble.feed.func<-spotFeedback.error.full} \cr
#' The function average model predictions, default is \code{spotConfig$seq.ensemble.average.func<-mean} \cr
#' The number of designs created for validation, which should be zero, i.e. \code{spotConfig$seq.ensemble.cut.num <- 0} 
#' The epsilon parameter, default is \code{spotConfig$seq.ensemble.epsilon <- 0} 
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
spotEnsembleMultiChoose <- function(rawB,mergedB,design,spotConfig,fit=NULL){	
	design <- spotInitializePredictor(design,"data.frame",spotConfig$alg.roi,NULL,"spotEnsembleMultiChoose",spotConfig$io.verbosity)
	########################################################################################
	# Initialize Settings (defaults)
	if(is.null(spotConfig$seq.ensemble.feed.func)){
		spotConfig$seq.ensemble.feed.func<-spotFeedback.error.full
	}	
	if(is.null(spotConfig$seq.ensemble.epsilon)){
		spotConfig$seq.ensemble.epsilon<-0 #zero means no randomization, one means completely random
	}
	#######################################################################################
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
		res <- spotEnsembleModelEvaluate(spotConfig,design,all=FALSE)	
	}
	#Check error
	if(is.na(err[1])){ #if no error from last step, calculate from this step (i.e. first step)
		spotConfig$seq.modelFit$fitList<-res$fits
		err=spotConfig$seq.ensemble.feed.func(spotConfig,mergedB,rawB)
	}
	#######################################################################################
	#choose a model
	if (runif(1)>spotConfig$seq.ensemble.epsilon){
		# Determine "best" model
		chosenMod=which(err==min(as.numeric(err)))#select mod with lowest error 
		chosenMod=chosenMod[1]; #if two or more models are the same, use first.
	}else{
		### choose random Model
		chosenMod <- sample(1:length(spotConfig$seq.ensemble.predictors),1);
	}   
	# set values to be returned inside spotConfig
	spotWriteLines(spotConfig$io.verbosity,3,"spotEnsembleMultiChoose finished");	
	spotConfig$seq.largeDesignY<-as.data.frame(res$Y[,1,chosenMod]);
	spotConfig$seq.prediction.prefModel<-c(spotConfig$seq.prediction.prefModel,chosenMod) #Remember which models were chosen
	spotConfig$seq.modelFit$err<-err
	spotConfig$seq.modelFit$fitList<-res$fits
	spotConfig$seq.modelFit$chosen<-chosenMod
	
	spotConfig
}

