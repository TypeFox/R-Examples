###################################################################################
#' Single Ensemble: Epsilon Greedy
#' 
#' Eps. Greedy - Select the model with the greatest success with probability 1-epsilon (greedy decision), 
#' select a random model with probability epsilon (explorative decision).
#' That means, the algorithm is completely greedy with epsilon=0, and completely explorative with epsilon=1\cr
#' The default reward function is \code{spotConfig$seq.ensemble.feed.func<-\link{spotFeedback.reward.bern}}.\cr
#' The value of epsilon \code{spotConfig$seq.greed.epsilon<-0.5}.
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
#' @references - Joannes Vermorel and Mehryar Mohri. 2005. Multi-armed bandit algorithms and empirical evaluation. 
#' In Proceedings of the 16th European conference on Machine Learning (ECML'05), Joao Gama, Rui Camacho, 
#' Pavel B. Brazdil, Alipio Mario Jorge, and Luis Torgo (Eds.). Springer-Verlag, Berlin, Heidelberg, 437-448.\cr
#' - M. Friese, M. Zaefferer, T. Bartz-Beielstein, O. Flasch, P. Koch, W. Konen, and
#' B. Naujoks. Ensemble based optimization and tuning algorithms. In F. Hoffmann
#' and E. Huellermeier, editors, Proceedings 21. Workshop Computational Intelligence,
#' p. 119-134. Universitaetsverlag Karlsruhe, 2011.
#'
#' @export
###################################################################################
spotEnsembleSingleEpsGreedy <- function(rawB,mergedB,design,spotConfig,fit=NULL){	
	spotWriteLines(spotConfig$io.verbosity,3,"spotEnsembleSingleEpsGreedy started");	
	models=spotConfig$seq.ensemble.predictors	
	########################################################
	# BUILD AND PREDICT
	########################################################		
	if(is.null(fit)){
		##############################DEFAULT HANDLING###################
		if(is.null(spotConfig$seq.ensemble.feed.func))
			spotConfig$seq.ensemble.feed.func<-spotFeedback.reward.bern;
		if(is.null(spotConfig$seq.greed.epsilon))
			spotConfig$seq.greed.epsilon<-0.5;
		#################################################################
		K=length(models)
		eps<-spotConfig$seq.greed.epsilon
		if(is.null(spotConfig$seq.greed.r)){#first initialization step		
			spotConfig$seq.greed.r<-rep(0,K)
			spotConfig$seq.greed.n<-rep(0,K)
		}	
		if(!is.null(spotConfig$seq.greed.imax)){ #calculate reward if a model was chosen in the last step
			imax<-spotConfig$seq.greed.imax
			reward<- spotConfig$seq.ensemble.feed.func(spotConfig,mergedB,rawB)
			if(reward>0){#success
				spotConfig$seq.greed.r[imax]=spotConfig$seq.greed.r[imax]+1
			}
			spotConfig$seq.greed.n[imax]<-spotConfig$seq.greed.n[imax]+1
		}
		meanReward<-spotConfig$seq.greed.r/spotConfig$seq.greed.n; 
		meanReward[which(is.na(meanReward))]<-0	#n is initialized with zero, so in the beginning mean values will be NaN	

		if(eps==1){ #if epsilon is one: completely random choice, (not excluding "best" model)
			imax<-sample(1:K,1)
		}else if(runif(1)<eps){
			imax<-sample((1:K)[-imax],1)#exploration decision
		}else{
			imax<-sample(which(meanReward==max(meanReward)),1)#greedy decision, sample if there is more than one maximum in reward vector
		}
		tmpConf <- append(spotConfig$seq.ensemble.settings[[imax]],spotConfig)
		conf <- eval(call(models[[imax]],rawB#evaluated chosen model
							, mergedB
							, design
							, tmpConf))
		spotConfig$seq.modelFit <- conf$seq.modelFit	
		spotConfig$seq.greed.imax<-imax#safe index, to associate reward with model
		spotConfig$seq.greed.history<-c(spotConfig$seq.greed.history,imax)
	}else{
		########################################################
		# PREDICT ONLY
		########################################################			
		jj<-spotConfig$seq.greed.imax
		tmpConf <- append(spotConfig$seq.ensemble.settings[[jj]],spotConfig)
		conf <- eval(call(models[[jj]], NULL
							, NULL
							, design
							, tmpConf
							, fit
							)); 
	}
	spotWriteLines(spotConfig$io.verbosity,3,"spotEnsembleSingleEpsGreedy finished");		
	spotConfig$seq.largeDesignY<-conf$seq.largeDesignY#usual output
	spotConfig
}