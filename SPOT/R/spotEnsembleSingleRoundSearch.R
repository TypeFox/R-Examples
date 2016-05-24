###################################################################################
#' Single Ensemble: RoundSearch
#' 
#' Round Search - This is a very simple approach to employing a set of models (ensemble).
#' Each of the chosen models is used in turn. That means:\cr
#' At time t, select model number \code{((t-1)mod k)+1}. \cr
#' In contrast to other single ensembles (i.e. \code{\link{spotEnsembleSingleSoftMax}}) 
#' this approach does not use any kind of reward function to rate the success of the models.
#' It has also no additional parameters to set, but will strongly depend on the ordering of the specified models.
#'
#' This is a "single ensemble", meaning that in every sequential step only one model in the ensemble is trained and evaluated.
#' The target is to actively "learn" which of the models are most suitable, based on their individual success.\cr
#' The models used are specified in the \code{spotConfig} list, for instance:\cr
#' \code{spotConfig$seq.ensemble.predictors = c(spotPredictRandomForest, spotPredictEarth, spotPredictForrester, spotPredictDace)}\cr
#' To specify the settings of each individual model, set:\cr
#' \code{seq.ensemble.settings = list(list(setting=1),list(setting=2),list(setting=3),list(setting=4))}\cr
#' Any parameters set in each of the corresponding lists (here: 4 individual lists) will overwrite settings in the main \code{spotConfig} list,
#' when the concerned model function is called.
#'
#' @param rawB matrix of raw x and y values
#' @param mergedB matrix of merged x and y values, does not have replicate entries
#' @param design design points to be evaluated by the meta model  
#' @param spotConfig the list of all parameters. Information about probability and success are stored in the sublist \code{spotConfig$seq.predictDual} 
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
####################################################################################
spotEnsembleSingleRoundSearch<- function(rawB,mergedB,design,spotConfig, fit=NULL){
	spotWriteLines(spotConfig$io.verbosity,3,"spotEnsembleSingleRoundSearch started");	
	########################################################
	# BUILD AND PREDICT
	########################################################			
	if(is.null(fit)){	
		### Choose next Predictor
		if (is.null(spotConfig$seq.predictDual$last)) spotConfig$seq.predictDual$last <- 0;
		chosenPredictor <- (spotConfig$seq.predictDual$last %% length(spotConfig$seq.ensemble.predictors)) + 1;
		
		### Get Prediction
		tmpConf <- append(spotConfig$seq.ensemble.settings[[chosenPredictor]],spotConfig)
		tmpRes <- eval(call(spotConfig$seq.ensemble.predictors[[chosenPredictor]]
							, rawB#evaluated chosen model
							, mergedB
							, design
							, tmpConf))
		
		### remember last Predictor
		spotConfig$seq.predictDual$last <- chosenPredictor;
		spotConfig$seq.modelFit<-tmpRes$seq.modelFit
	} else{ 
		########################################################
		# PREDICT ONLY
		########################################################			
		chosenPredictor<-spotConfig$seq.predictDual$last
		tmpConf <- append(spotConfig$seq.ensemble.settings[[chosenPredictor]],spotConfig)
		tmpRes <- eval(call(spotConfig$seq.ensemble.predictors[[chosenPredictor]],NULL,NULL, design, tmpConf ,fit));
	}
	spotWriteLines(spotConfig$io.verbosity,3,"spotEnsembleSingleRoundSearch finished");		
	spotConfig$seq.largeDesignY <- tmpRes$seq.largeDesignY;	#return predicted value	
	spotConfig
}
