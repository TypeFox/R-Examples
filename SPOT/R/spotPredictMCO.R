###################################################################################
#' Meta Model Interface: Multi Criteria Modelling
#' 
#' This interface function is supposed to be used for Multi Criteria Problems. It can be employed when the user
#' wants to specify different models for each of the objectives, instead of modelling all the objectives with the same technique.
#' The user has therefore to specify a list of configurations, where the different models and their settings are specified.
#'
#' @param rawB unmerged data
#' @param mergedB merged data
#' @param design new design points which should be predicted
#' @param spotConfig global list of all options, needed to provide data for calling functions. The most important elements in this list are here:\cr
#'	\code{spotConfig$mco.configs} list of model configurations, e.g. \code{=list(list(seq.predictionModel.func="spotPredictForrester",seq.forr.lambda.upval=-5),list(seq.predictionModel.func="spotPredictForrester",seq.forr.lambda.upval=-1))}\cr
#'	In this example, two Kriging models are specified for each of two objectives, but with different settings for the lower boundary of lambda. Else, different models could be specified, e.g., \code{=list(list(seq.predictionModel.func="spotPredictForrester"),list(seq.predictionModel.func="spotPredictLm"))}
##' @param fit if an existing model fit is supplied, the model will not be build based on 
#'				data, but only evaluated with the model fit (on the design data). To build the model, 
#'				this parameter has to be NULL. If it is not NULL the parameters mergedB and rawB will not be 
#'				used at all in the function.
#'
#' @return returns the list \code{spotConfig} with two new entries:\cr
#' 	spotConfig$seq.modelFit fit of the model used with the predictor functions \cr
#'	spotConfig$seq.largeDesignY the y values of the design, evaluated with the fit
#' @export
###################################################################################
spotPredictMCO  <- function(rawB,mergedB,design,spotConfig,fit=NULL){	
	#design <- spotInitializePredictor(design,"data.frame",spotConfig$alg.roi,NULL,"spotPredictMCO",spotConfig$io.verbosity)	
	spotWriteLines(spotConfig$io.verbosity,3,"spotPredictMCO started");
	xNames <- row.names(spotConfig$alg.roi);
	yNames <- spotConfig$alg.resultColumn	
	########################################################
	# BUILD
	########################################################
	if(is.null(fit)){#build and predict	
		spotConfig$seq.modelFit.y<-rawB[yNames]
		fit=list()
		resy=NULL		
		resvar=NULL		
		for (i in 1:length(yNames)){
			tmpConf <- append(spotConfig$mco.configs[[i]],spotConfig)
			tmpConf <- tmpConf[unique(names(tmpConf))]
			tmpConf$alg.resultColumn <- yNames[i]
			tmpConf$io.verbosity <- spotConfig$io.verbosity
			tmp<- eval(call(tmpConf$seq.predictionModel.func
							, rawB[c(yNames[i],xNames)]
							, mergedB[c(yNames[i],xNames)]
							, design
							, tmpConf))
			spotConfig$mco.configs[[i]] = tmp
			resy <- cbind(resy,unlist(tmp$seq.largeDesignY))
			if(spotConfig$seq.model.variance) resvar <- cbind(resvar,unlist(tmp$seq.largeDesignVar))
			fit[[i]] <- tmp$seq.modelFit
		}			
	}else{ #predict only
		resy=NULL
		resvar=NULL
		for (i in 1:length(yNames)){
			tmpConf <- append(spotConfig$mco.configs[[i]],spotConfig)
			tmpConf <- tmpConf[unique(names(tmpConf))]
			tmpConf$alg.resultColumn <- yNames[i]	
			tmpConf$io.verbosity <- spotConfig$io.verbosity			
			tmp<- eval(call(tmpConf$seq.predictionModel.func
							, NULL
							, NULL
							, design
							, tmpConf
							, spotConfig$seq.modelFit[[i]]))
			resy <- cbind(resy,unlist(tmp$seq.largeDesignY))
			if(spotConfig$seq.model.variance) resvar <- cbind(resvar,unlist(tmp$seq.largeDesignVar))
		}			
	}	
	if(is.function(spotConfig$seq.infill) & !is.null(design)){# do EI 
		y=spotConfig$alg.currentResult[spotConfig$alg.resultColumn] 
		resy= spotConfig$seq.infill(resy,resvar,y,spotConfig$mco.refPoint)
	}	
	########################################################
	# OUTPUT
	########################################################	
	spotWriteLines(spotConfig$io.verbosity,3,"spotPredictMCO finished")
	spotConfig$seq.largeDesignY=as.data.frame(resy)	
	spotConfig$seq.modelFit<-fit;
	spotConfig
}

