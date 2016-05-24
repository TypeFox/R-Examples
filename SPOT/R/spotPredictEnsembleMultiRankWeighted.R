###################################################################################
#' spotEnsembleMultiRankWeighted
#' 
#' Alternating Recommendation of Candidate Points \cr
#' This "multi ensemble" uses the same ranking approach as \code{\link{spotEnsembleMultiRank}}. Additionally, the ranks
#' are weighted. To determine an appropriate weight, the models feedback (e.g. error measure) is employed.\cr
#' Relevant configuration parameters of this ensemble are: \cr
#' The function to calculate feedback (e.g. the error) of the models, default is \code{spotConfig$seq.ensemble.feed.func<-spotFeedback.error.full} \cr
#' The number of designs cut off for leave-one-out validation should always be zero for this ensemble, i.e \code{spotConfig$seq.ensemble.cut.num <- 0} 
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
#' @param fit should always be NULL. Evaluation of existing fits is not implemented for this ensemble.
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
spotEnsembleMultiRankWeighted <- function(rawB,mergedB,design,spotConfig,fit=NULL){
	design <- spotInitializePredictor(design,"data.frame",spotConfig$alg.roi,NULL,"spotEnsembleMultiRankWeighted",spotConfig$io.verbosity)
	### rate last predictions
	#MSPEs <- calcPredictionErrorsForAllModels(spotConfig, mergedB, design);
	
	### Fit all Models and get Predictions		
	#allPredictions <- fitAllModels(spotConfig$ensemble.predictors, numberOfFittings, cut.type, cut.size, rawB, mergedB, design, spotConfig);
	
	
	# Initialize
	if(is.null(spotConfig$seq.ensemble.feed.func)){
		spotConfig$seq.ensemble.feed.func<-spotFeedback.error.full #set default function
	}
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
		stop("Plotting or optimization of existing model ensembles of this type is not implemented.") #prob 
	}	

	### do a voting
	votedPrediction <- spotEnsembleSummedRankWeighted(res, err, spotConfig$seq.design.size);
	
	spotWriteLines(spotConfig$io.verbosity,3,"spotEnsembleMultiRankWeighted finished");	
	### return votings (as largeDesignY)	
	spotConfig$seq.largeDesignY<-as.data.frame(votedPrediction);
	spotConfig$seq.modelFit$err<-err
	spotConfig$seq.modelFit$fitList<-res$fits
	spotConfig
}

