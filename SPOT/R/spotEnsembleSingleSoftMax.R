###################################################################################
#' Single Ensemble: SoftMax
#' 
#' SoftMax Algorithm - Select a model based on a probability proportional to its success (i.e. reward) \cr
#' The default reward function is \code{spotConfig$seq.ensemble.feed.func<-\link{spotFeedback.reward.bern}}. \cr
#' Two parameters have to be set. The parameter alpha is used as a multiplier when a models reward is updated:\cr
#' The default value is \code{spotConfig$seq.softmax.alpha <- 0.5}. \cr
#' The parameter tau is used to control exploration/exploitation trade-off. SoftMax is completely greedy with tau=0, and completely random with tau going to infinity:\cr
#' The default value is \code{spotConfig$seq.softmax.tau <- 1}. 
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
#' @references - Sutton, R. S.; Barto, A. G.: Reinforcement Learning: An Introduction (Adaptive Computation
#' and Machine Learning). MIT Press. URL http://www.cse.iitm.ac.in/~cs670/
#' book/the-book.html. 1998.\cr
#'- Joannes Vermorel and Mehryar Mohri. 2005. Multi-armed bandit algorithms and empirical evaluation. 
#' In Proceedings of the 16th European conference on Machine Learning (ECML'05), Joao Gama, Rui Camacho, 
#' Pavel B. Brazdil, Alipio Mario Jorge, and Luis Torgo (Eds.). Springer-Verlag, Berlin, Heidelberg, 437-448.\cr
#' - M. Friese, M. Zaefferer, T. Bartz-Beielstein, O. Flasch, P. Koch, W. Konen, and
#' B. Naujoks. Ensemble based optimization and tuning algorithms. In F. Hoffmann
#' and E. Huellermeier, editors, Proceedings 21. Workshop Computational Intelligence,
#' p. 119-134. Universitaetsverlag Karlsruhe, 2011.
#'
#' @export
###################################################################################
spotEnsembleSingleSoftMax<- function(rawB,mergedB,design,spotConfig,fit=NULL){
	spotWriteLines(spotConfig$io.verbosity,3,"spotEnsembleSingleSoftMax started");	
	models=spotConfig$seq.ensemble.predictors
	########################################################
	# BUILD AND PREDICT
	########################################################		
	if(is.null(fit)){
		##############################DEFAULT HANDLING###################
		if(is.null(spotConfig$seq.ensemble.feed.func))
			spotConfig$seq.ensemble.feed.func<-spotFeedback.reward.bern;
		#################################################################	
		K=length(models)
		alpha<-spotConfig$seq.softmax.alpha
		tau<-spotConfig$seq.softmax.tau
		if(is.null(alpha)){alpha<-0.5}
		if(is.null(tau)){tau<-1}
		if(is.null(spotConfig$seq.softmax.q)){
			#first initialization step
			spotConfig$seq.softmax.q<-rep(0,K)
			spotConfig$seq.softmax.p<-rep(1/K,K)
		}
		if(!is.null(spotConfig$seq.softmax.imax)){ #calculate reward if a model was chosen in the last step
			imax<-spotConfig$seq.softmax.imax
			reward<-spotConfig$seq.ensemble.feed.func(spotConfig,mergedB,rawB) #TODO
			spotConfig$seq.softmax.q[imax]<-spotConfig$seq.softmax.q[imax]+alpha*(reward-spotConfig$seq.softmax.q[imax])
			for(j in 1:K){
				spotConfig$seq.softmax.p[j]=exp(spotConfig$seq.softmax.q[j]/tau)/sum(exp(spotConfig$seq.softmax.q/tau))
			}
		}
		imax<-sum(cumsum(spotConfig$seq.softmax.p)<=runif(1))+1	
		tmpConf <- append(spotConfig$seq.ensemble.settings[[imax]],spotConfig)
		conf <- eval(call(models[[imax]],rawB#evaluated chosen model
							, mergedB
							, design
							, tmpConf))
		spotConfig$seq.modelFit<-conf$seq.modelFit#safe the model
		spotConfig$seq.softmax.imax<-imax#safe index, to associate reward with model
		spotConfig$seq.softmax.history<-c(spotConfig$seq.softmax.history,imax)
	}else{ #Else do only prediction on preexisting information
		########################################################
		# PREDICT ONLY
		########################################################			
		jj<-spotConfig$seq.softmax.imax
		tmpConf <- append(spotConfig$seq.ensemble.settings[[jj]],spotConfig)
		conf <- eval(call(models[[jj]], NULL
							, NULL
							, design
							, tmpConf
							, fit
							)); 
	}	
	spotWriteLines(spotConfig$io.verbosity,3,"spotEnsembleSingleSoftMax finished");		
	spotConfig$seq.largeDesignY<-conf$seq.largeDesignY#usual output	
	spotConfig
}