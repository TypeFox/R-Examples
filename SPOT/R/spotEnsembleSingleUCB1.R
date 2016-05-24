###################################################################################
#' Single Ensemble: UCB1
#' 
#' UCB1 - Upper Confidence Bound 1 \cr
#' Use the model that maximizes \code{xi + sqrt( 2 * ln(n) / ni  )}, where code{xi} is the average reward of the model \code{i}, 
#' \code{n} is the number of sequential steps (i.e. number of times any model was used), and \code{ni} is the number of times the model \code{i} was chosen.
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
#' @references - Peter Auer, Nicolo Cesa-Bianchi, and Paul Fischer. Finite-time analysis of the multiarmed
#' bandit problem. Machine Learning, 47((2-3)), 2002.\cr
#' - Volodymyr Kuleshov and Doina Precup. Algorithms for the multi-armed bandit problem.
#' Submitted, 2010. URL: http://www.cs.mcgill.ca/~vkules/bandits.pdf\cr
#' - M. Friese, M. Zaefferer, T. Bartz-Beielstein, O. Flasch, P. Koch, W. Konen, and
#' B. Naujoks. Ensemble based optimization and tuning algorithms. In F. Hoffmann
#' and E. Huellermeier, editors, Proceedings 21. Workshop Computational Intelligence,
#' p. 119-134. Universitaetsverlag Karlsruhe, 2011.
#'
#' @export
###################################################################################
spotEnsembleSingleUCB1<- function(rawB,mergedB,design,spotConfig,fit=NULL){	
	spotWriteLines(spotConfig$io.verbosity,3,"spotEnsembleSingleUCB1 started");
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
		if(is.null(spotConfig$seq.ucb1.n)){
			#first initialization step
			spotConfig$seq.ucb1.mj<-rep(1,K)
			spotConfig$seq.ucb1.nj<-rep(1,K)
			spotConfig$seq.ucb1.n<-1
		}
		if(!is.null(spotConfig$seq.ucb1.imax)){ #calculate reward if a model was chosen in the last step
			imax<-spotConfig$seq.ucb1.imax
			reward<-spotConfig$seq.ensemble.feed.func(spotConfig,mergedB,rawB)
			
			spotConfig$seq.ucb1.nj[imax]<-spotConfig$seq.ucb1.nj[imax]+1
			spotConfig$seq.ucb1.n<-spotConfig$seq.ucb1.n+1
			if(reward>0){spotConfig$seq.ucb1.mj[imax]<-spotConfig$seq.ucb1.mj[imax]+1}#success
		}
		x<-rep(0,K)
		for(i in 1:K){ #take a sample from beta distribution for each model, to select the model 
			x[i]<-(spotConfig$seq.ucb1.mj[i]/spotConfig$seq.ucb1.nj[i])+sqrt((2*log(spotConfig$seq.ucb1.n))/(spotConfig$seq.ucb1.nj[i]))
		}
		imax<- which.max(x) #argmax(x)	
		tmpConf <- append(spotConfig$seq.ensemble.settings[[imax]],spotConfig)
		conf <- eval(call(models[[imax]],rawB#evaluated chosen model
							, mergedB
							, design
							, tmpConf))
		spotConfig$seq.modelFit<-conf$seq.modelFit#safe the model
		spotConfig$seq.ucb1.imax<-imax#safe index, to associate reward with model
		spotConfig$seq.ucb1.history<-c(spotConfig$seq.ucb1.history,imax)
	}else{ #Else do only prediction on preexisting information
		########################################################
		# PREDICT ONLY
		########################################################		
		jj<-spotConfig$seq.ucb1.imax
		tmpConf <- append(spotConfig$seq.ensemble.settings[[jj]],spotConfig)
		conf <- eval(call(models[[jj]], NULL
							, NULL
							, design
							, tmpConf
							, fit
							)); 
	}	
	spotWriteLines(spotConfig$io.verbosity,3,"spotEnsembleSingleUCB1 finished");	
	spotConfig$seq.largeDesignY<-conf$seq.largeDesignY#usual output	
	spotConfig
}
