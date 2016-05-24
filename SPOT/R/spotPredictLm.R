###################################################################################
#' Meta Model Interface: Linear Model
#'
#' A linear prediction model, which will use higher order interactions if data is sufficient.
#' Can be used both for single and multi objective SPOT.
#' 
#' This function implements a linear model for prediction. Depending on the numbers of variables 
#' either no interactions,  interaction  between the variables may be used or a full quadratic model is
#' provided.  
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
spotPredictLm <- function(rawB,mergedB,design,spotConfig,fit=NULL){	
	design <- spotInitializePredictor(design,"data.frame",spotConfig$alg.roi,"rsm","spotPredictLm",spotConfig$io.verbosity)
	
	########################################################
	# BUILD
	########################################################	

	if(is.null(fit)){
		nParam <- nrow(spotConfig$alg.roi)	
		pNames <- row.names(spotConfig$alg.roi)	
		makeNNames <- function(n) { Map(function(i) paste("x", i, sep = ""), 1:n) }
		codeNames <- makeNNames(nParam)
		codeNames[spotConfig$alg.roi$type == "FACTOR"] <- pNames[spotConfig$alg.roi$type == "FACTOR"] #Factors. not. coded.
		yNames <- spotConfig$alg.resultColumn
		X <- rawB[pNames] #MZ: Bugfix for 1 dimensional optimization
		if(length(yNames)>1){ #for multi objective modeling
			fit<-list()
			res=list()
			for(ii in 1:length(yNames)){
				Y<-rawB[yNames[ii]]					
				df1 <- data.frame(cbind(Y,X))			
				fmla <- spotCodeFormulas(X,pNames,spotConfig$alg.roi)				
				rsmDf <- coded.data(df1, formulas = fmla)
				#check feasibility of all points: use only those IN the ROI
				rsmDf <- rsmDf[apply(rsmDf[,-1], 1, function(x) all(x >= -1 & x <= 1)), ]
				## Check if there are sufficient observations for linear, quadratic etc. models
				nExp <- nrow(X)
				## Linear model without interactions:
				nRequired1 <- 1 + nParam
				## Linear model with interactions:
				nRequired2 <- 1 + nParam * (nParam + 1) / 2 
				## Full quadratic model:
				nRequired3 <- 1 + nParam + nParam * (nParam + 1) / 2
				## Create the most powerful model possible given the current number of parameters...
				paramString <- paste(codeNames,collapse=",") # the string "x1,x2,...,xnParams"
				if (nExp >= nRequired1 && nExp < nRequired2) {
					spotWriteLines(spotConfig$io.verbosity,3,"spotPredictLm: First order (FO) effects estimated by rsm.");
					rsmFormula <- as.formula(sprintf(paste(colnames(Y),"~ FO(%s)")
												, paramString))
					fit1 <- rsm(formula = rsmFormula, data = rsmDf)
				}else if (nExp >= nRequired2 && nExp < nRequired3) {
					spotWriteLines(spotConfig$io.verbosity,3,"spotPredictLm: First order (FO) with two-way interactions (TWI) effects estimated by rsm.");
					rsmFormula <- as.formula(sprintf(paste(colnames(Y),"~ FO(%s) + TWI(%s)")
												, paramString, paramString))
					fit1 <- rsm(formula = rsmFormula, data = rsmDf)
				}else if (nExp >= nRequired3) {
					spotWriteLines(spotConfig$io.verbosity,3,"spotPredictLm: Second order (SO) effects estimated by rsm.");
					rsmFormula <- as.formula(sprintf(paste(colnames(Y),"~ FO(%s) + TWI(%s) + PQ(%s)") ## rsmFormula <- as.formula(sprintf("Y ~ SO(%s)", paramString)) generates warnings, so we use the long form
														, paramString, paramString, paramString))
					fit1 <- rsm(formula = rsmFormula, data = rsmDf)
				}
				fit[[ii]]<-fit1
			}
		}else{#for single objective modeling
			Y<-rawB[yNames]
			df1 <- data.frame(cbind(Y,X))			
			fmla <- spotCodeFormulas(X,pNames,spotConfig$alg.roi)		
			rsmDf <- coded.data(df1, formulas = fmla)
			#check feasibility of all points: use only those IN the ROI
			rsmDf <- rsmDf[apply(rsmDf[,-1], 1, function(x) all(x >= -1 & x <= 1)), ]
			###
			## Check if there are sufficiently many data for linear, quadratic etc. models
			nExp <- nrow(X)
			## Linear model without interactions:
			nRequired1 <- 1 + nParam
			## Linear model with interactions:
			nRequired2 <- 1 + nParam * (nParam + 1) / 2 
			## Full quadratic model:
			nRequired3 <- 1 + nParam + nParam * (nParam + 1) / 2
			## Create the most powerful model possible given the current number of parameters...
			paramString <-  paste(codeNames,collapse=",") # the string "x1,x2,...,xnParams"
			if (nExp >= nRequired1 && nExp < nRequired2) {
				spotWriteLines(spotConfig$io.verbosity,3,"spotPredictLm: First order (FO) effects estimated by rsm.");
				rsmFormula <- as.formula(sprintf(paste(colnames(Y),"~ FO(%s)")
												, paramString))
				fit <- rsm(formula = rsmFormula, data = rsmDf)
			}else if (nExp >= nRequired2 && nExp < nRequired3) {
				spotWriteLines(spotConfig$io.verbosity,3,"spotPredictLm: First order (FO) with two-way interactions (TWI) effects estimated by rsm.");
				rsmFormula <- as.formula(sprintf(paste(colnames(Y),"~ FO(%s) + TWI(%s)")
												, paramString, paramString))
				fit <- rsm(formula = rsmFormula, data = rsmDf)
			}else if (nExp >= nRequired3) {
				spotWriteLines(spotConfig$io.verbosity,3,"spotPredictLm: Second order (SO) effects estimated by rsm.");
				rsmFormula <- as.formula(sprintf(paste(colnames(Y),"~ FO(%s) + TWI(%s) + PQ(%s)")## rsmFormula <- as.formula(sprintf("Y ~ SO(%s)", paramString)) generates warnings, so we use the long form
														, paramString, paramString, paramString))
				fit <- rsm(formula = rsmFormula, data = rsmDf)
			}
		}
	}
	########################################################
	# PREDICT
	########################################################	
	if(!is.null(design)){
		if(length(spotConfig$alg.resultColumn)>1){
			res=list()
			for (i in 1:length(fit)){
				codings <- fit[[i]]$coding
				design <- val2code(design, codings)
				res[[i]]<-predict(fit[[i]],design)
			}	
		}
		else{
			codings <- fit$coding
			design <- val2code(design, codings)
			res <- predict(fit,design)
		}	
	}else{res <- NULL}
	########################################################
	# OUTPUT
	########################################################	
	spotWriteLines(spotConfig$io.verbosity,3,"spotPredictLm finished")
  	spotConfig$seq.modelFit<-fit
	spotConfig$seq.largeDesignY<-as.data.frame(res)
	spotConfig
}
