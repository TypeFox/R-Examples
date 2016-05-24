###################################################################################
#' Multi criteria optimization of predicted surrogate models
#'  
#' Uses by default the number of design points expected as population size for multi criteria optimization of the
#' models build in the current sequential step. Executed after building the prediction models.
#' 
#' @param startPoint initial information, not yet used
#' @param spotConfig list of all options, needed to provide data for calling functions.
#' 
#' This function uses the parameter \code{spotConfig$seq.modelFit}. This is supposed to be a list of fits (i.e. one fit for each objective), 
#' that can be evaluated by calling the original model interface with that fit list.
#' The parameter \code{spotConfig$seq.predictionOpt.method} will be used to choose the optimization method to be used to find the minimum of the fitted model:\cr
#' "nsga2" \cr
#' "sms-emoa"\cr
#' 
#' @return returns the list \code{spotConfig} with two new entries:\cr
#' 	spotConfig$optDesign are the parameters of the new Pareto optimal design points \cr
#'	spotConfig$optDesignY are the associated values of the objective functions (e.g. meta model values)
#' @seealso  \code{\link{spotModelOptim}}solves the same task for single objective optimization (i.e. just one surrogate model)\cr
#' See \code{\link{spotSmsEmoa}} for the used SMS-EMOA implementation
#' @export
###################################################################################
spotModelParetoOptim <- function(startPoint,spotConfig){
	#startPoint=as.numeric(startPoint[1,])
	spotWriteLines(spotConfig$io.verbosity,3,"spotModelParetoOptim started")
	if(is.null(spotConfig$seq.predictionOpt.method)) spotConfig$seq.predictionOpt.method = "sms-emoa"
	if(is.null(spotConfig$seq.predictionOpt.budget)) spotConfig$seq.predictionOpt.budget = 200
	if(is.null(spotConfig$seq.predictionOpt.psize)) spotConfig$seq.predictionOpt.psize = 20
	if(is.null(spotConfig$seq.predictionOpt.psize)) spotConfig$seq.predictionOpt.sbx.n = 15
	if(is.null(spotConfig$seq.predictionOpt.psize)) spotConfig$seq.predictionOpt.sbx.p = 0.7
	if(is.null(spotConfig$seq.predictionOpt.psize)) spotConfig$seq.predictionOpt.pm.n = 25
	if(is.null(spotConfig$seq.predictionOpt.psize)) spotConfig$seq.predictionOpt.pm.p = 0.3
	if(is.null(spotConfig$seq.predictionOpt.restarts)) spotConfig$seq.predictionOpt.restarts = FALSE
	fit<-spotConfig$seq.modelFit;	
	tempVerbose<-spotConfig$io.verbosity
	#building the local objective function
	spotGetFitness <- function(x){
		spotConfig$io.verbosity=0;
		as.numeric(eval(call(spotConfig$seq.predictionModel.func
								, NULL 
								, NULL
								, as.data.frame(x)
								, spotConfig
                , fit #external fit is used, model is only evaluated not build, therefore the NULLS are no prob
								))$seq.largeDesignY)
	}
	lowROI<-as.numeric(spotConfig$alg.aroi[[1]])
	upROI<-as.numeric(spotConfig$alg.aroi[[2]])
	subMethod=spotConfig$seq.predictionOpt.method
	budget=spotConfig$seq.predictionOpt.budget
	psize=spotConfig$seq.predictionOpt.psize
	if(is.null(psize))max(2,ceiling(0.1*budget))  #limit population size to number of points to be evaluated!
	
	
	r1 <- spotOptimizationInterfaceMco(startPoint,spotGetFitness,gr=NULL,lower=lowROI,upper=upROI,method=subMethod,control=list(
			fevals = spotConfig$seq.predictionOpt.budget
			,popsize = spotConfig$seq.predictionOpt.psize
			,restarts = spotConfig$seq.predictionOpt.restarts
			,sbx.n=spotConfig$seq.predictionOpt.sbx.n 
			,sbx.p=spotConfig$seq.predictionOpt.sbx.p
			,pm.n=spotConfig$seq.predictionOpt.pm.n 
			,pm.p=spotConfig$seq.predictionOpt.pm.p
			),ref=rep(NA,length(spotConfig$alg.resultColumn))) #todo: REFpoint
	newDesignPrediction <- r1$value
	newDesign <- r1$par


	spotPrint(spotConfig$io.verbosity,3,newDesign)
	spotWriteLines(spotConfig$io.verbosity,3,"spotModelParetoOptim finished")
	
	spotConfig$optDesign<-newDesign
	spotConfig$optDesignY<-newDesignPrediction
	spotConfig$io.verbosity<-tempVerbose
	spotConfig
}

###################################################################################
#' DEPRECATED
#'  
#' see \code{\link{spotModelParetoOptim}}
#'
#' @param startPoint see \code{\link{spotModelParetoOptim}}
#' @param spotConfig see \code{\link{spotModelParetoOptim}}
#' 
#' @export
#' @keywords internal
###################################################################################
spotParetoOptMulti <- function(startPoint,spotConfig){
	warning("spotParetoOptMulti is a deprecated function, please use spotModelParetoOptim")
	spotModelParetoOptim(startPoint,spotConfig)
}