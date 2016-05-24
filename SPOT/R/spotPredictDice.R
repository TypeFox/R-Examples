###################################################################################
#' Meta Model Interface: Dice Kriging  
#'
#' Kriging meta model based on the DiceKriging package.
#' It usually provides good prediction performance, but is rather unstable. 
#' 
#' @param rawB matrix of raw x and y values
#' @param mergedB matrix of merged x and y values, does not have replicate entries
#' @param design design points to be evaluated by the meta model 
#' @param spotConfig the list of all parameters is given.
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
####################################################################################
spotPredictDice <- function(rawB,mergedB,design,spotConfig,fit=NULL){	
	design <- spotInitializePredictor(design,"data.frame",spotConfig$alg.roi,"DiceKriging","spotPredictDice",spotConfig$io.verbosity)	
	########################################################
	# BUILD
	########################################################	
	if(is.null(fit)){
		xNames <- row.names(spotConfig$alg.roi)
		yNames <- spotConfig$alg.resultColumn
		x <- rawB[xNames]
		if(length(yNames)==1){
			y <- data.frame(y=rawB[[yNames]])
			if(spotConfig$io.verbosity>2){fit<-DiceKriging::km(design=x,response=y,nugget.estim=TRUE)}
			else{fit<-DiceKriging::km(design=x,response=y,control=list(trace=FALSE),nugget.estim=TRUE)}
			res <- predict(fit,design,"UK")$mean 
		}
		else{#Distinction for multi criteria 
			y <- data.frame(y=rawB[yNames])
			fit=list()
			res=list()
			for (i in 1:length(yNames)){
				if(spotConfig$io.verbosity>2){fit[[i]]<-DiceKriging::km(design=x,response=y[,i],nugget.estim=TRUE)}
				else{fit[[i]]<-DiceKriging::km(design=x,response=y[,i],control=list(trace=FALSE),nugget.estim=TRUE)}
				res[[i]]<-predict(fit[[i]],design,"UK")$mean 
			}			
		}		
	}else{ #use existing model!
		fit<-fit
	}
	########################################################
	# PREDICT
	########################################################	
	if(!is.null(design)){ 	
		if(length(spotConfig$alg.resultColumn)>1){
			res=list()
			for (i in 1:length(fit)){
				res[[i]]<-predict(fit[[i]],design,"UK")$mean 
			}	
		}
		else{res <- predict(fit,design,"UK")$mean}			
	}else{res <- NULL}			
	########################################################
	# OUTPUT
	########################################################	
	spotWriteLines(spotConfig$io.verbosity,3,"spotPredictDice finished")
	spotConfig$seq.modelFit<-fit
	spotConfig$seq.largeDesignY<-as.data.frame(res)
	spotConfig
}
