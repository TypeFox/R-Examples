###################################################################################
#' Meta Model Interface: Maximum Likelihood Estimation for Gaussian Processes, Kriging 
#'
#' Kriging model based on mlegp package.
#' This function uses two settings, which are stored in the spotConfig parameter:\cr
#' \code{spotConfig$seq.mlegp.constantMean} Use constant mean (mu) in mlegp (=1) or linear model (=0); 1 by default\cr
#' \code{spotConfig$seq.mlegp.min.nugget} minimum value of nugget term; 0 by default\cr\cr
#' If those settings are not in spotConfig their mentioned defaults will be used. \cr\cr
#' If the numeric value of \code{spotConfig$mlegp.reduce} is smaller than the observations in \code{mergedB}, 
#' \code{spotConfig$mlegp.reduce} will specify how many samples should be drawn without replacement from mergedB. 
#' This can prevent explosion of time consumption in this function.
#' Mlegp can be used both for single and multi objective SPOT.
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
spotPredictMlegp <- function(rawB,mergedB,design,spotConfig,fit=NULL){
	design <- spotInitializePredictor(design,"data.frame",spotConfig$alg.roi,"mlegp","spotPredictMlegp",spotConfig$io.verbosity)	
	########################################################
	# BUILD
	########################################################	
	if(is.null(fit)){
		if(is.null(spotConfig$seq.mlegp.constantMean)) spotConfig$seq.mlegp.constantMean <- 1 #default handling for user options  
		if(is.null(spotConfig$seq.mlegp.min.nugget)) spotConfig$seq.mlegp.min.nugget <- 0.0       	
		xNames <- row.names(spotConfig$alg.roi)
		yNames <- spotConfig$alg.resultColumn		
		x <- unname(as.matrix(rawB[xNames]))
		y <- unname(as.matrix(rawB[yNames]))
#####################################################################
		constantMean <- spotConfig$seq.mlegp.constantMean
		if (constantMean != 1) {
			ones <- rep(1, dim(x)[1])
			dx <- cbind(ones, x)
			t <- try(solve(t(dx) %*% dx), TRUE)
			if (class(t) == "try-error") {
				constantMean <-1
			}		
		}
#######################################################################################################		
#		START ERROR Handling and Model Building	
#######################################################################################################		
		fit<-mlegp::mlegp(X=x,Z=y, verbose = spotConfig$io.verbosity
                      , constantMean = constantMean
                      , min.nugget = spotConfig$seq.mlegp.min.nugget)
		if(is.null(fit) || any(sapply(fit,is.null))){ #mlegp will create a NULL fit sometimes, seems to be fixable by using non zero min.nugget
			#browser()
			fit<-mlegp::mlegp(X=x,Z=y, verbose = spotConfig$io.verbosity
                      , constantMean = constantMean
                      , min.nugget = 1)
		}
		#TODO:Problem: adapt this second problem for mco
#######################################################################################################	
#		END ERROR Handling and Model Building		
#######################################################################################################				
	}
	########################################################
	# PREDICT
	########################################################	
	if(!is.null(design)){ 	
		if(length(spotConfig$alg.resultColumn)==1){
			res<-predict(fit,design)
		}
		else{#Distinction for multi criteria spot 
			res=list()
			for (i in 1:length(spotConfig$alg.resultColumn)){
				res[[i]]<-predict(fit[[i]],design)
			}			
		}
	}else{res <- NULL}		
	#
	spotWriteLines(spotConfig$io.verbosity,2,"spotPredictMlegp finished successfully")
	spotConfig$seq.modelFit<-fit
	spotConfig$seq.largeDesignY<-as.data.frame(res)
	spotConfig
}
