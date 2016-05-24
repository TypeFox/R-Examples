###################################################################################
#' Single Ensemble: POKER
#'
#' POKER - Price Of Knowledge and Estimated Reward \cr
#' This approach uses the POKER algorithm to balance cost of knowledge acquisition against the
#' estimated reward. The single model with the highest rating is chosen. Success (i.e. reward)
#' is rated based on the improvement yielded by the model. If no new best point is found, 
#' reward should be zero.\cr
#' The default reward function is \code{spotConfig$seq.ensemble.feed.func<-\link{spotFeedback.reward.norm}}.\cr
#' The default horizon of the POKER algorithm (H) is \code{spotConfig$seq.poker.H<-10}. 
#'
#' This is a "single ensemble", meaning that in every sequential step only one model in the ensemble is trained and evaluated.
#' The target is to actively "learn" which of the models are most suitable, based on their individual success.\cr
#' The models used are specified in the \code{spotConfig} list, for instance:\cr
#' \code{spotConfig$seq.ensemble.predictors = c(spotPredictRandomForest, spotPredictEarth, spotPredictForrester, spotPredictDace)}\cr
#' To specify the settings of each individual models, set:\cr
#' \code{seq.ensemble.settings = list(list(setting=1),list(setting=2),list(setting=3),list(setting=4))}\cr
#' Any parameters set in each of the corresponding lists (here: 4 individual lists) will overwrite settings in the main \code{spotConfig} list,
#' when the concerned model function is called.
#'
#' @note Implementation might still be faulty, not suggested for serious experiments.
#' This is work in progress.
#'
#' @param rawB unmerged data
#' @param mergedB merged data
#' @param design new design points which should be predicted
#' @param spotConfig global list of all options, needed to provide data for calling functions
#' @param fit if an existing model fit is supplied, the model will not be build based on 
#'				data, but only evaluated with the model fit (on the design data). To build the model, 
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
spotEnsembleSinglePOKER <- function(rawB,mergedB,design,spotConfig,fit=NULL){
	spotWriteLines(spotConfig$io.verbosity,3,"spotEnsembleSinglePOKER started");
	models=spotConfig$seq.ensemble.predictors	
	########################################################
	# BUILD AND PREDICT
	########################################################		
	if(is.null(fit)){
		##############################DEFAULT HANDLING###################
		if(is.null(spotConfig$seq.ensemble.feed.func))
			spotConfig$seq.ensemble.feed.func<-spotFeedback.reward.norm;	
		if(is.null(spotConfig$seq.poker.H)) #horizon of the POKER algorithm
			spotConfig$seq.poker.H<-10;
		#################################################################		
		K=length(models)
		#if(is.null(spotConfig$seq.poker.T)){spotConfig$seq.poker.T<-10}
		H<-spotConfig$seq.poker.H #exploration parameter	
		#################################################################		
		#First calculate rewards if possible	
		if(!is.null(spotConfig$seq.poker.imax)){ #calculate reward if a model was chosen in the last step
			#reward is calculated by difference most recent "best" values
			reward<-spotConfig$seq.ensemble.feed.func(spotConfig,mergedB,rawB)
			spotConfig$seq.poker.r[spotConfig$seq.poker.imax]<-spotConfig$seq.poker.r[spotConfig$seq.poker.imax]+reward
			spotConfig$seq.poker.r2[spotConfig$seq.poker.imax]<-spotConfig$seq.poker.r2[spotConfig$seq.poker.imax]+reward^2
		}
		
		if(is.null(spotConfig$seq.poker.n)){
		#first initialization step
			spotConfig$seq.poker.n<-rep(1,K)
			spotConfig$seq.poker.r<-rep(1,K)
			spotConfig$seq.poker.r2<-rep(1,K)
			spotConfig$seq.poker.t<-1
			imax<-sample(1:K,1)
		}
		else if(spotConfig$seq.poker.t==2){
		#second initialization step, still completely random
			temp<-1:K
			temp<-temp[-spotConfig$seq.poker.imax]
			if(length(temp)==1){imax<-temp}
			else{imax<-sample(temp,1)}
		}
		else{
			q<- length(which(spotConfig$seq.poker.r>0))
			i0<-which.max(spotConfig$seq.poker.r/spotConfig$seq.poker.n)
			temp <- max(round(sqrt(q)),1)
			i1<-which(temp==order(spotConfig$seq.poker.r))
			deltamu<-((spotConfig$seq.poker.r[i0]/spotConfig$seq.poker.n[i0])-(spotConfig$seq.poker.r[i1]/spotConfig$seq.poker.n[i1]))/round(sqrt(q))
			pmax<-(-Inf)#initialize pmax and imax
			imax<-NULL
			for(i in 1:K){
				if(spotConfig$seq.poker.n[i]>0){
					n<-spotConfig$seq.poker.n[i]
					mu<-spotConfig$seq.poker.r[i]/spotConfig$seq.poker.n[i]
				}else{
					n<-1 #better solution? else there will be divided by zero for p<-mu+deltamu*.......			
					k<-which(spotConfig$seq.poker.n>0)
					mu<-sum(spotConfig$seq.poker.r[k]/spotConfig$seq.poker.n[k])/sum(spotConfig$seq.poker.n[k])			
				}
				if(spotConfig$seq.poker.n[i]>1){
					sigma<-sqrt((spotConfig$seq.poker.r2[i]/spotConfig$seq.poker.n[i])-((spotConfig$seq.poker.r[i]/spotConfig$seq.poker.n[i])^2))
				}else{
					k<-which(spotConfig$seq.poker.n>1)
					sigma<-sum(sqrt((spotConfig$seq.poker.r2[k]/spotConfig$seq.poker.n[k])-((spotConfig$seq.poker.r[k]/spotConfig$seq.poker.n[k])^2)))/sum(spotConfig$seq.poker.n[k])		
				}	
				#hier H durch T-t ersetzen?
				p<-mu+deltamu*H*(1-pnorm(i0+deltamu,mu,sigma/sqrt(n)))
				if(p>pmax){
					pmax<-p
					imax<-i
				}
			}	
		}
		tmpConf <- append(spotConfig$seq.ensemble.settings[[imax]],spotConfig)
		conf <- eval(call(models[[imax]],rawB#evaluated chosen model
							, mergedB
							, design
							, tmpConf))
		spotConfig$seq.modelFit<-conf$seq.modelFit						
		spotConfig$seq.poker.imax<-imax#safe index, to associate reward with model
		spotConfig$seq.poker.n[imax]<-spotConfig$seq.poker.n[imax]+1#increase number of "lever pulls" for chosen model
		spotConfig$seq.poker.t<-spotConfig$seq.poker.t+1#increment internal step counter						
	}else{ #Else do only prediction on preexisting information
		########################################################
		# PREDICT ONLY
		########################################################	
		imax<-spotConfig$seq.poker.imax
		tmpConf <- append(spotConfig$seq.ensemble.settings[[imax]],spotConfig)
		conf <- eval(call(models[[imax]], NULL
							, NULL
							, design
							, tmpConf
							, fit
							)); 
	}	
	spotWriteLines(spotConfig$io.verbosity,3,"spotEnsembleSinglePOKER finished");	
	spotConfig$seq.largeDesignY<-conf$seq.largeDesignY	#usual output
	spotConfig
}
