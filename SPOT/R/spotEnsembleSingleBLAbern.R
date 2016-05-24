###################################################################################
#' Single Ensemble: BLAbern
#' 
#' BLA Bernoulli - Bayesian Learning Automaton for Bernoulli Distributed Feedback[Braadland_Norheim]
#' An advantage of this is, that the only relevant parameters are the initial values of the beta distribution.
#' A disadvantage is the stationarity of the algorithm. The ensemble problem in SPOT is of dynamic nature.\cr
#' The default reward function is \code{spotConfig$seq.ensemble.feed.func<-\link{spotFeedback.reward.bern}}.
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
#' @references - O.-C. Granmo. A Bayesian Learning Automaton for Solving Two-Armed Bernoulli 
#' Bandit Problems. Machine Learning and Applications, ICMLA '08. p. 23-30. 2008.\cr
#' - M. Friese, M. Zaefferer, T. Bartz-Beielstein, O. Flasch, P. Koch, W. Konen, and
#' B. Naujoks. Ensemble based optimization and tuning algorithms. In F. Hoffmann
#' and E. Huellermeier, editors, Proceedings 21. Workshop Computational Intelligence,
#' p. 119-134. Universitaetsverlag Karlsruhe. 2011.
#'
#' @export
###################################################################################
spotEnsembleSingleBLAbern<- function(rawB,mergedB,design,spotConfig,fit=NULL){
	spotWriteLines(spotConfig$io.verbosity,3,"spotEnsembleSingleBLAbern started")
	models=spotConfig$seq.ensemble.predictors
	########################################################
	# BUILD AND PREDICT
	########################################################	
	if(is.null(fit)){
		K=length(models)
		if(is.null(spotConfig$seq.bla.a)){
			#first initialization step
			spotConfig$seq.bla.a<-rep(1,K)
			spotConfig$seq.bla.b<-rep(1,K)
		}
		if(!is.null(spotConfig$seq.bla.imax)){ #calculate reward if a model was chosen in the last step
			imax<-spotConfig$seq.bla.imax
			reward<- spotFeedback.reward.bern(spotConfig,mergedB,rawB)
			if(reward<=0){spotConfig$seq.bla.b[imax]=spotConfig$seq.bla.b[imax]+1}#failure
			else{spotConfig$seq.bla.a[imax]=spotConfig$seq.bla.a[imax]+1}		#success
		}
		x<-rep(0,K)
		for(i in 1:K){ #take a sample from beta distribution for each model, to select the model 
			x[i]<-rbeta(1,spotConfig$seq.bla.a[i],spotConfig$seq.bla.b[i])
		}
		imax<- which.max(x) #argmax(x)	
		tmpConf <- append(spotConfig$seq.ensemble.settings[[imax]],spotConfig)
		conf <- eval(call(models[[imax]],rawB#evaluated chosen model
							, mergedB
							, design
							, tmpConf))
		spotConfig$seq.modelFit<-conf$seq.modelFit#safe the model
		spotConfig$seq.bla.imax<-imax#safe index, to associate reward with model
		spotConfig$seq.bla.history<-c(spotConfig$seq.bla.history,imax)
	}else{ #Else do only prediction on pre-existing information
		########################################################
		# PREDICT ONLY
		########################################################			
		jj<-spotConfig$seq.bla.imax
		tmpConf <- append(spotConfig$seq.ensemble.settings[[jj]],spotConfig)
		conf <- eval(call(models[[jj]], NULL
							, NULL
							, design
							, tmpConf
							, fit
							)); 
	}	
	spotWriteLines(spotConfig$io.verbosity,3,"spotEnsembleSingleBLAbern finished")	
	spotConfig$seq.largeDesignY<-conf$seq.largeDesignY#usual output	
	spotConfig
}
