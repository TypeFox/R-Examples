###################################################################################
#' Meta Model Interface: Random Forest
#'  
#' A prediction model interface based on randomForest package, using a random forest
#' for regression.
#' Can be used both for single and multi objective SPOT.
#'
#' This is a model that can incorporate parameters which are marked as FACTORS (i.e. categorical parameters) in the region of interest, see \code{\link{spotROI}}.
#' Please note that the design used to train the RF model should contain all levels of the factor variable. FACTORS are not ordered, and therefore are impossible to extrapolate on.
#' If new data is given in the \code{design} variable which contains unseen FACTOR levels, please note that this will probably create NA values in the prediction. 
#' NA values might yield errors in your SPOT run, ending it prematurely. It is therefore recommended to build a initial design which contains at least one example
#' of each FACTOR level.
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
#' @return returns the list \code{spotConfig} with two new entries:\cr
#' 	spotConfig$seq.modelFit fit of the model used with predict() \cr
#'	spotConfig$seq.largeDesignY the y values of the design, evaluated with the fit
#'
#' @seealso \code{\link{SPOT}}
#' @export
###################################################################################
spotPredictRandomForest <- function(rawB,mergedB,design,spotConfig,fit=NULL){	
	design <- spotInitializePredictor(design,"data.frame",spotConfig$alg.roi,"randomForest","spotPredictRandomForest",spotConfig$io.verbosity)	
	xNames <- row.names(spotConfig$alg.roi)
	########################################################
	# BUILD
	########################################################
	if(is.null(fit)){#IF NO EXTERNAL FIT: BUILD MODEL AND EVALUATE		
		yNames <- spotConfig$alg.resultColumn
		x <- rawB[xNames]
		#force observation-locations into factors 
		tmp <- spotForceFactorX(x,spotConfig$alg.roi$type,xNames)
		x <- tmp$x #parameters are now factors if specified in roi
		spotConfig$alg.factor.levels <- tmp$levels #needed for later conversion of prediction locations		
		if(length(yNames)==1){
			spotConfig$seq.modelFit.min<-min(rawB[yNames])
			y <- rawB[[yNames]]
			fit <- randomForest(x, y)		
		}else{#Distinction for multi criteria spot 
			y <- rawB[yNames]
			spotConfig$seq.modelFit.y<-rawB[yNames]
			fit=list()
			for (i in 1:length(yNames)){
				fit[[i]]<-randomForest(x,y[,i])
			}			
		}
	}
	########################################################
	# PREDICT
	########################################################
	if(!is.null(design)){ 
		pred.all=spotConfig$seq.model.variance
		nmodel <- length(spotConfig$alg.resultColumn)
		design <- spotForceFactorDesign(design,spotConfig$alg.roi$type,xNames,spotConfig$alg.factor.levels)
		if(nmodel>1){ #do multi criteria prediction
			resy=matrix(0,nrow(design),nmodel)
			resvar=matrix(NA,nrow(design),nmodel)
			y=list()
			for (i in 1:length(fit)){ #predict			
				if(pred.all){		
					res <- predict(fit[[i]],design,predict.all=TRUE)
					resy[,i]= res$aggregate
					resvar[,i]= apply(res$individual,1,sd)#use prediction from each tree to calc variance
				}else{
					resy[,i]= predict(fit[[i]],design)
				}		
				y[[i]]= spotConfig$seq.modelFit.y[,i]
			}
			if(is.function(spotConfig$seq.infill)){# do EI 
				resy= spotConfig$seq.infill(resy,resvar,y,spotConfig$mco.refPoint)
			}
		}else{ #do single criteria prediction
			resvar=matrix(NA,nrow(design),1)
			if(pred.all){		
				res <- predict(fit,design,predict.all=TRUE)
				resy <- res$aggregate
				resvar <- apply(res$individual,1,sd) #use prediction from each tree to calc variance
			}else{
				resy <- predict(fit,design)
			}		
			if(is.function(spotConfig$seq.infill)){ # do EI		
				resy= spotConfig$seq.infill(resy,resvar,spotConfig$seq.modelFit.min)
			}
		}	
	}else{
		resy <- NULL
		resvar <- NULL
	}		
	########################################################
	# OUTPUT
	########################################################	
	spotWriteLines(spotConfig$io.verbosity,3,"spotPredictRandomForest finished")
	spotConfig$seq.largeDesignY=as.data.frame(resy)	
	spotConfig$seq.largeDesignVar=as.data.frame(resvar)	
	spotConfig$seq.modelFit<-fit
	spotConfig
}

