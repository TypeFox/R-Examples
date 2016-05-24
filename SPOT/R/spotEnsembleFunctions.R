#############################################################################
#############################################################################
#THIS ARE ALL FEEDBACK FUNCTIONS FOR ENSEMBLES. reward, or error... etc.
#############################################################################
#############################################################################
###################################################################################################
#' Ensemble Feedback Functions
#' 
#' This describes the feedback functions which can be employed in the \code{spotEnsembleSingle*} and \code{spotEnsembleMulti*} functions.
#'
#' If \code{spotFeedback.error.full} or \code{spotFeedback.error.last} are used, the user can choose an error function of type \code{f(a,b)} where \code{a}
#' and \code{b} are vectors whose differences are the deviations between model and real function. By default this function is chosen with
#' \code{spotConfig$seq.ensemble.error.func <- spotSelectionRmse}. See \code{\link{spotSelectionCriteria}} for some other options.
#' Note that ".full" indicates that the modelling error is calculated based on all known design points, whereas ".last" indicates that the modelling error is only calculated based 
#' on the last evaluated design points, which can be considered unseen data for the model. Therefore, ".full" might be biased, whereas ".last" might be more random due to the small
#' number of samples. \code{spotFeedback.error.full} and \code{spotFeedback.error.last} return scalar values for each model, i.e. a vector.\cr\cr
#' \code{spotFeedback.y} simply returns predictions from previous steps. It therefore returns a vector for each model, i.e. a matrix of num-model columns, and npoints rows.\cr\cr
#' \code{spotFeedback.sd.interSubModelsFull} when cut off sub-designs are available for validation, this function can be used to calculate the standard deviation at each point. It therefore returns a vector for each model, i.e. a matrix of num-model columns, and npoints rows.\cr\cr
#' \code{spotFeedback.reward*} reward based on success of the model used in last sequential step. Scalar value returned. 
#'
#' @name spotFeedback
#' @usage spotFeedback.y(spotConfig,mergedB,rawB)
#' spotFeedback.sd.interSubModelsFull(spotConfig,mergedB,rawB)
#' spotFeedback.deviation(spotConfig,mergedB,rawB)
#' spotFeedback.error.full(spotConfig,mergedB,rawB)
#' spotFeedback.error.last(spotConfig,mergedB,rawB)
#' spotFeedback.error.order(spotConfig,mergedB,rawB)
#' spotFeedback.error.combo(spotConfig,mergedB,rawB)
#' spotFeedback.reward.bern(spotConfig,mergedB,rawB)
#' spotFeedback.reward.norm(spotConfig,mergedB,rawB)
#' @param spotConfig parameter list, as created by the calling functions
#' @param mergedB merged list of design points as evaluated on the target function of SPOT
#' @param rawB raw list of design points as evaluated on the target function of SPOT
#'
#' @return \code{spotFeedback.y, spotFeedback.sd.interSubModelsFull and spotFeedback.deviation} return a matrix with as many columns as models in the ensemble and as many rows as points in mergedB. \cr
#'  \code{spotFeedback.reward.bern and spotFeedback.reward.norm} return single values, which indicate the Bernoulli or normal distributed reward for the model used in the last sequential step. \cr
#' all other functions return a vector with an error measure for each model in the ensemble.
#'
#' @export spotFeedback.y spotFeedback.sd.interSubModelsFull spotFeedback.deviation spotFeedback.error.full spotFeedback.error.last spotFeedback.error.order spotFeedback.error.combo spotFeedback.error.combo spotFeedback.reward.bern spotFeedback.reward.norm
#' @aliases spotFeedback.y spotFeedback.sd.interSubModelsFull spotFeedback.deviation spotFeedback.error.full spotFeedback.error.last spotFeedback.error.order spotFeedback.error.combo spotFeedback.error.combo spotFeedback.reward.bern spotFeedback.reward.norm
###################################################################################################
# just return the predictions from last step
spotFeedback.y <- function(spotConfig,mergedB,rawB){
  models=spotConfig$seq.ensemble.predictors;
  nmodels=length(models);
    
  # prepare empty ResultContainer
  results=as.data.frame(array(0, c(nrow(mergedB[1]),nmodels)));
  
  for (i in 1:nmodels){ ### For every predictor
    spotConfig1<-spotConfig;
    spotConfig1$io.verbosity=0;
    spotConfig1 <- eval(call(models[[i]],NULL, NULL, as.data.frame(mergedB[-1])
                               ,spotConfig1
                               ,spotConfig$seq.modelFit$fitList[[i]][[1]]))
    
    results[,i] <- spotConfig1$seq.largeDesignY ### redo last Fit
  }
  results
  
}

# for each modeltype, calc the sd per point over all submodels 
spotFeedback.sd.interSubModelsFull <- function(spotConfig,mergedB,rawB){
  models=spotConfig$seq.ensemble.predictors;
  nmodels=length(models);
  
  # prepare empty ResultContainer
  results=as.data.frame(array(0, c(nrow(mergedB[1]),nmodels)));
  
  for (i in 1:nmodels){ ### For every predictor
    spotConfig1<-spotConfig;
    spotConfig1$io.verbosity=0;
    countSubmodels <- length(spotConfig$seq.modelFit$fitList[[1]]);
    subResults=as.data.frame(array(0, c(nrow(mergedB[1]),countSubmodels)));    
    
    for (sub in 1:countSubmodels){
      subModel <- eval(call(models[[i]],NULL, NULL, as.data.frame(mergedB[-1])
                     ,spotConfig1
                     ,spotConfig$seq.modelFit$fitList[[i]][[sub]]));
      subResults[,sub] <- subModel$seq.largeDesignY; ### redo last Fit      
    }    
    results[,i] <- apply(subResults, 1, sd)
  }  
  results
}

# returns matrix, with errors at each point for each model  (matrix is number of models times number of points )
spotFeedback.deviation <- function(spotConfig,mergedB,rawB){
  m=list();
  models=spotConfig$seq.ensemble.predictors;
  nmodels=length(models);
    
  # prepare empty ResultContainer
  results=array(0, c(nrow(mergedB[1]),nmodels));                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
  
	for (i in 1:nmodels){ ### For every predictor
    spotConfig1<-spotConfig;
  	spotConfig1$io.verbosity=0;
		spotConfig1 <- eval(call(models[[i]],NULL, NULL, as.data.frame(mergedB[-1])
                                ,spotConfig1
								,spotConfig$seq.modelFit$fitList[[i]][[1]]))
    m[[i]]=spotConfig1$seq.largeDesignY ### redo last Fit
    results[,i] <- abs((m[i] - mergedB[1])[,1]) ### calc absolute prediction Error
	}
  results
}

#returns vector with overall error measure for each model (number of models is length of returned vector)
spotFeedback.error.full <- function(spotConfig,mergedB,rawB){ 
	#m=list();
	err=array();
	models=spotConfig$seq.ensemble.predictors
	nmodels=length(models)
	if(is.null(spotConfig$seq.ensemble.error.func))spotConfig$seq.ensemble.error.func <- spotSelectionRmse
	for (i in 1:nmodels){
		spotConfig1<-spotConfig;
		spotConfig1$io.verbosity=0;
		spotConfig1 <- eval(call(models[[i]],NULL, NULL, as.data.frame(mergedB[-1])
                                ,spotConfig1
								,spotConfig$seq.modelFit$fitList[[i]][[1]]));
		#m[[i]]=spotConfig1$seq.largeDesignY;
		err[i] <- spotConfig$seq.ensemble.error.func(as.numeric(unlist(spotConfig1$seq.largeDesignY)),  as.numeric(unlist(mergedB[1])))
	}
	err
}	


#returns vector with overall error measure for each model (number of models is length of returned vector)
spotFeedback.error.last <- function(spotConfig,mergedB,rawB){ 
	#m=list();
	err=array();
	models=spotConfig$seq.ensemble.predictors
	nmodels=length(models)
	mergedData <- spotPrepareData(spotConfig)
	if(is.null(spotConfig$seq.ensemble.error.func))spotConfig$seq.ensemble.error.func <- spotSelectionRmse	
	for (i in 1:nmodels){
		spotConfig1<-spotConfig;
		spotConfig1$io.verbosity=0;
		spotConfig1 <- eval(call(models[[i]],NULL, NULL, as.data.frame(mergedData$x[mergedData$STEP==max(mergedData$STEP),])
                                ,spotConfig1
								,spotConfig$seq.modelFit$fitList[[i]][[1]]));
		#m[[i]]=spotConfig1$seq.largeDesignY;		
		err[i] <- spotConfig$seq.ensemble.error.func(as.numeric(unlist(spotConfig1$seq.largeDesignY)),  as.numeric(unlist(mergedData$mergedY[mergedData$STEP==max(mergedData$STEP)])))
	}
	err
}	

#Das ist: Ordnungsfehler des Modelles, ueber alle Punkte
spotFeedback.error.order  <-  function(spotConfig,mergedB,rawB){  
	m=list()
	err=array()
	models=spotConfig$seq.ensemble.predictors
	nmodels=length(models)
	for (i in 1:nmodels){
		spotConfig1<-spotConfig
		spotConfig1$io.verbosity=0
		spotConfig1 <- eval(call(models[[i]],NULL, NULL, as.data.frame(mergedB[-1])
                                ,spotConfig1
								,spotConfig$seq.modelFit$fitList[[i]][[1]]))
		m[[i]]=spotConfig1$seq.largeDesignY
		err[i]=sum(abs(order(m[[i]])-order(mergedB[1])))/length(mergedB[1])
	}
	err
}

#Das ist: Kombination aus full und order error
spotFeedback.error.combo <-  function(spotConfig,mergedB,rawB){ 
	err<-spotFeedback.error.order(spotConfig,mergedB,rawB)
	err2<-spotFeedback.error.full(spotConfig,mergedB,rawB)
	err<-((err/max(err))+(err2/max(err2)))/2
	chosenMod=which(err==min(as.numeric(err)))				
	chosenMod=chosenMod[1]; 
}

#Das ist: calculate bernoulli distrib. feedback for single approach ensembles, success 1 or failure 0
spotFeedback.reward.bern <-  function(spotConfig,mergedB,rawB){ 
	spotConfig$alg.currentBest[nrow(spotConfig$alg.currentBest),]$CONFIG-spotConfig$alg.currentBest[nrow(spotConfig$alg.currentBest)-1,]$CONFIG
}

#Das ist: calculate normal distrib. feedback for single approach ensembles
spotFeedback.reward.norm <-  function(spotConfig,mergedB,rawB){ 
	if(spotConfig$alg.currentBest[nrow(spotConfig$alg.currentBest),]$CONFIG>spotConfig$alg.currentBest[nrow(spotConfig$alg.currentBest)-1,]$CONFIG){
		reward<-abs(spotConfig$alg.currentBest[nrow(spotConfig$alg.currentBest),1]-spotConfig$alg.currentBest[nrow(spotConfig$alg.currentBest)-1,1])
	}else{reward=0}	
	reward
}

#############################################################################
#############################################################################
#next: Function(s) for model building/eval
#############################################################################
#############################################################################

###################################################################################
#' spotEnsembleModelBuilding
#' 
#' This function is used in the multi ensembles (e.g. \code{\link{spotEnsembleMultiAlternate}}) to
#' build and evaluate the models in the ensemble. 
#'
#' @param rawB unmerged data, passed to the model function (i.e. \code{spotPredict*})
#' @param mergedB merged data, passed to the model function (i.e. \code{spotPredict*}
#' @param design new design points which should be predicted
#' @param spotConfig global list of all options
#' 
#' @return returns a list, containing the element \code{Y} which is a list of results from each mode, and the element \code{fits} which is a list of model fits.
#' @keywords internal
###################################################################################
spotEnsembleModelBuilding <- function( rawB, mergedB, design, spotConfig){	
  
	models=spotConfig$seq.ensemble.predictors #Models in ensemble
	configs=spotConfig$seq.ensemble.settings #list of configuration for each model in ensemble
	if(is.null(spotConfig$seq.ensemble.cut.num)) spotConfig$seq.ensemble.cut.num = 0 #default number of subsample designs. 0 means no subdesigns!
	results <- list(); #initialize returned list
	results$Y <- array(0,c(nrow(design),spotConfig$seq.ensemble.cut.num+1,length(models))) #3d array consisting of as many 2d matrixes as models exist
	if(spotConfig$seq.ensemble.cut.num>0){	
		subs <- spotEnsembleSubdesigns(spotConfig$alg.currentResult,row.names(spotConfig$alg.roi),spotConfig$alg.resultColumn,
										spotConfig$seq.ensemble.cut.size,spotConfig$seq.ensemble.cut.type,spotConfig$seq.ensemble.cut.num,
										spotConfig$seq.transformation.func,spotConfig$seq.merge.func) #genereate designs on subsamples of data
	}
	for (i in 1:length(models)){	
		### predict on full design:
		tmpConf=append(configs[[i]],spotConfig)
		tmpRes <- eval(call(models[[i]],rawB, mergedB, design, tmpConf))
		results$Y[,1,i]<-unlist(tmpRes$seq.largeDesignY) #save Y values
		results$fits[[i]] <- list() #save fit
		results$fits[[i]][[1]] <- tmpRes$seq.modelFit
		### predict on subdesigns:
		if(spotConfig$seq.ensemble.cut.num>0){
			for (n in 1:length(subs$rawB)){			
				tmpRes <- eval(call(models[[i]],subs$rawB[[n]], subs$mergB[[n]], design, tmpConf))
				results$Y[,n+1,i] <- unlist(tmpRes$seq.largeDesignY)
				results$fits[[i]][[n+1]] <- tmpRes$seq.modelFit
			}	
		}			
	}	
	results
}

###################################################################################
#' spotEnsembleSubdesigns
#' 
#' Cut off smaller subdesigns for ensembles
#'
#' @param res result data
#' @param xnames parameter names
#' @param ynames objective names
#' @param cut.size size of the sub designs
#' @param cut.type type of cut-off, "perc" is percentage, otherwise cut.size is an absolute number of points cut off
#' @param cut.num number of subdesigns created
#' 
#' @return returns a list, containing the element \code{rawB} with raw design and \code{mergB} with merged design
#' @keywords internal
###################################################################################
spotEnsembleSubdesigns <- function(res,xnames,ynames,cut.size=10,cut.type="perc",cut.num=1,merge,transform){	
	if(is.null(cut.type)) cut.type = "perc"
	if(is.null(cut.size)) cut.size = 10.0

	design.size <- nrow(res)
	 
	### cut defines the number of points that will be dismissed for each subdesign
	cut <- if (cut.type=="perc") round(design.size*cut.size/100) else cut.size
	if (cut<1)cut <- 1; ### at least one point has to be dismissed

	designs <- NULL
	designsMerged <- NULL
	
	for (i in 1:cut.num){ #remember: taking out one point in rawB also affects mergedB! hence both are reconstructed here.
		cutThis <- sample(1:design.size, cut, replace=FALSE);
		des<- res[-cutThis,]                   ### reduced subdesign
		z <- split(des[,ynames], des$CONFIG);
		fs <- function(xs) { transform(merge(xs)) }
		mergedY <- sapply(z,fs)
		if (length(xnames)==1){ 
			mergedX <- as.data.frame(sapply(split(des[,xnames], des$CONFIG),mean))
			names(mergedX)<-xnames
		}else{
			mergedX <- as.data.frame(t(sapply(split(des[,xnames], des$CONFIG),colMeans)))
		}
		y <- des[,ynames]
 		B <- data.frame(res[order(y,decreasing=FALSE),])
		mergedB <- data.frame(cbind(mergedY,as.matrix(mergedX))[order(mergedY,decreasing=FALSE),])
		designs <- if (is.null(designs)) list(B) else c(designs, list(B));		
		designsMerged <- if (is.null(designsMerged)) list(mergedB) else c(designsMerged, list(mergedB))
	}	
	list(rawB=designs,mergB=designsMerged)
}

###################################################################################
#' spotEnsembleModelEvaluate
#' 
#' This function is used in the multi ensembles (e.g. \code{\link{spotEnsembleMultiAlternate}}) to
#' evaluate the models in the ensemble. 
#'
#' @param spotConfig global list of all options
#' @param design new design points which should be predicted
#' @param all if false, only the chosen model is evaluated, else all models are evaluated
#' 
#' @return returns a list, containing the element \code{Y} which is a list of results from each mode, and the element \code{fits} which is a list of model fits.
#' @keywords internal
###################################################################################
spotEnsembleModelEvaluate <- function(spotConfig,design,all){	
	res<-list()	
	models=spotConfig$seq.ensemble.predictors #Models in ensemble
	configs=spotConfig$seq.ensemble.settings #list of configuration for each model in ensemble
	design<-data.frame(design)	
	if(ncol(design)!=nrow(spotConfig$alg.roi)){ #ugly, but necessary because R is strange about converting one dimensional vectors to dataframes (confusing rows/columns)
		design<-t(design)			
	}

	res$Y <- array(0,c(nrow(design),1,length(models)))
		
	
	if(all==FALSE){ #only the chosen model is evaluated
		jj<-spotConfig$seq.modelFit$chosen		
		conf <- eval(call(models[[jj]], NULL
							, NULL
							, design
							, append(configs[[jj]],spotConfig) #use additional settings 
							, spotConfig$seq.modelFit$fitList[[jj]][[1]])); 
		res$Y[,1,jj]=unlist(conf$seq.largeDesignY);		
	}
	else{ #evaluate all models 
		for (jj in length(models):1){	
			conf <- eval(call(models[[jj]], NULL
								, NULL
								, design
								, append(configs[[jj]],spotConfig) #use additional settings 
								, spotConfig$seq.modelFit$fitList[[jj]][[1]])); 
			res$Y[,1,jj]=unlist(conf$seq.largeDesignY);				
		}		
	}
	res$fits=spotConfig$seq.modelFit$fitList #fits remain unchanged
	res
}

###################################################################################
#' spotEnsembleStatsMin
#' 
#' Get minimum of ensemble results.
#'
#' @param x result structure: a list with element \code{Y} for predicted values and \code{fits} for model fit list
#' 
#' @return returns a vector of results, with one element of the vector for each fit
#' @keywords internal
###################################################################################
spotEnsembleStatsMin <- function(x){
	results=array(1,length(x$fits))
	for (i in 1:length(x)){
		results[i]= min(x$Y[[1]])
	}	
	results
}

###################################################################################
#' spotEnsembleStatsSd
#' 
#' Get standard deviation for each point.
#'
#' @param x result structure: a list with element \code{Y} for predicted values and \code{fits} for model fit list
#' 
#' @return returns a matrix of results, with one column for each fit and one row for each pont
#' @keywords internal
###################################################################################
spotEnsembleStatsSd <- function(x){ 
	results=array(0, c(nrow(x$Y[,,1]),length(x$fits)))
	for (i in 1:length(x$fits)){
		results[,i]=  apply(t(x$Y[,,i]),2,sd)
	}	
	results
}

###################################################################################
#' spotEnsembleSummedRank
#' 
#' Get the voting result for each point. Summed over all models.
#'
#' @param x result structure: a list with element \code{Y} for predicted values and \code{fits} for model fit list
#' @param size size of design to be modeled
#' 
#' @return the summed voted order
#' @keywords internal
###################################################################################
spotEnsembleSummedRank <- function(x, size){
	countModels <- length(x$fits);	
	voteOrder <- rep(0,size);
	
	for (i in 1:countModels){		
		predictions <- x$Y[,1,i]#x[[i]][,1];
		ranked <- rank(predictions);
		voteOrder <- voteOrder + (exp(ranked));
	}		
	voteOrder
}

###################################################################################
#' spotEnsembleMinimum 
#' 
#' Minn over all points.
#'
#' @param x result structure: a list with element \code{Y} for predicted values and \code{fits} for model fit list
#' 
#' @return the summed 
#' @keywords internal
###################################################################################
spotEnsembleMinimum <- function(x){
  apply(as.data.frame(x$Y), 1, min)
}

####################################################################################
#' spotEnsembleSummedRankWeighted
#' 
#' Get the voting result for each point. Summed over all models.
#'
#' @param x result structure: a list with element \code{Y} for predicted values and \code{fits} for model fit list
#' @param error mean squared prediction errors
#' @param size size of design to be modeled
#' 
#' @return the summed voted order, weighted by error
#' @keywords internal
###################################################################################
spotEnsembleSummedRankWeighted <- function(x, error, size){
	weight <- (error/max(error+1))+1; ## scaling to ]1,2[
	countModels <- length(x$fits);	
	voteOrder <- rep(0,size);
	
	for (i in 1:countModels){		
		predictions <- x$Y[,1,i]#x[[i]][,1];
		ranked <- rank(predictions)/weight[i];
		voteOrder <- voteOrder + (exp(ranked));
	}		
	voteOrder		
}
