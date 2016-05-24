###################################################################################
#' Single Ensemble: BLAnorm
#' 
#' BLA normal - Bayesian Learning Automaton for Normally Distributed Feedback[Braadland_Norheim]
#' An advantage of this approach is, that the only relevant parameters are the initial values which here remain at defaults.
#' A disadvantage is the stationarity of the algorithm. The ensemble problem in SPOT is of dynamic nature.\cr
#' The default reward function is \code{spotConfig$seq.ensemble.feed.func<-\link{spotFeedback.reward.norm}}.
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
spotEnsembleSingleBLAnorm<- function(rawB,mergedB,design,spotConfig,fit=NULL){	
	spotWriteLines(spotConfig$io.verbosity,3,"spotEnsembleSingleBLAnorm started");	
	models=spotConfig$seq.ensemble.predictors
	########################################################
	# BUILD AND PREDICT
	########################################################		
	if(is.null(fit)){
		##############################DEFAULT HANDLING###################
		if(is.null(spotConfig$seq.ensemble.feed.func))
			spotConfig$seq.ensemble.feed.func<-spotFeedback.reward.norm;
		#################################################################
		K=length(models)
		if(is.null(spotConfig$seq.bla.mu)){
			#first initialization step
			spotConfig$seq.bla.mu<-rep(0,K)
			spotConfig$seq.bla.v<-rep(1,K)#degrees of freedom		
			spotConfig$seq.bla.k<-rep(1,K)#number of observations
			spotConfig$seq.bla.vsd<-rep(1,K)
		}	
		if(!is.null(spotConfig$seq.bla.imax)){ #calculate reward if a model was chosen in the last step
			#reward is calculated by difference most recent "best" values. If 
			i<-spotConfig$seq.bla.imax
			reward <- spotConfig$seq.ensemble.feed.func(spotConfig,mergedB,rawB)
			spotConfig$seq.bla.vsd[i]<-spotConfig$seq.bla.vsd[i]+(spotConfig$seq.bla.k[i]/(spotConfig$seq.bla.k[i]+1))*(reward-spotConfig$seq.bla.mu[i])^2
			spotConfig$seq.bla.mu[i]<-(spotConfig$seq.bla.k[i]/(spotConfig$seq.bla.k[i]+1))*spotConfig$seq.bla.mu[i]+(reward/(spotConfig$seq.bla.k[i]+1))
			spotConfig$seq.bla.v[i]<-spotConfig$seq.bla.v[i]+1	
			spotConfig$seq.bla.k[i]<-spotConfig$seq.bla.k[i]+1
		}	
		x<-rep(0,K)
		for(i in 1:K){ #take a sample from normalized inverse chi square distribution for each model, to select the model 
			S2<-spotInvChisq(1,spotConfig$seq.bla.v[i], spotConfig$seq.bla.vsd[i]/spotConfig$seq.bla.v[i])
			x[i]<-rnorm(1,spotConfig$seq.bla.v[i], S2/spotConfig$seq.bla.k[i])
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
	}else{ #Else do only prediction on preexisting information
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
	spotWriteLines(spotConfig$io.verbosity,3,"spotEnsembleSingleBLAnorm finished")
	spotConfig$seq.largeDesignY<-conf$seq.largeDesignY#usual output	
	spotConfig
}

########################################################		
#' Inverse Chisquared distribution
#'
#' @param n number of observations
#' @param df degrees of freedom
#' @param scale scale
#' @return Inverse Chisquare
#' @keywords internal
########################################################		
spotInvChisq<-function (n, df, scale){
    (df * scale)/rchisq(n, df = df)
}