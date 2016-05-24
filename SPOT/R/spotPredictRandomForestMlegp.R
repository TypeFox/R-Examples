###################################################################################
#' Meta Model Interface: Random Forest combined with Mlegp
#'  
#' A very simple ensemble which uses results from a Gaussian process model (mlegp) and a random forest.
#' 
#' @param rawB unmerged data
#' @param mergedB merged data
#' @param largeDesign new design points which should be predicted
#' @param spotConfig global list of all options, needed to provide data for calling functions
#' @param fit if an existing model fit is supplied, the model will not be build based on 
#'				data, but only evaluated with the model fit (on the largeDesign data). To build the model, 
#'				this parameter has to be NULL. If it is not NULL the parameters mergedB and rawB will not be 
#'				used at all in the function.
#'
#' @return returns the list \code{spotConfig} with two new entries:\cr
#' 	spotConfig$seq.modelFit fit of the model used with predict() \cr
#'	spotConfig$seq.largeDesignY the y values of the large design, evaluated with the fit
#'
#' @seealso \code{\link{SPOT}}
#' @export
###################################################################################
spotPredictRandomForestMlegp <- function(rawB,mergedB,largeDesign,spotConfig,fit=NULL){	     
	spotWriteLines(spotConfig$io.verbosity,3,"spotPredictRandomForestMlegp started");

	if(is.null(fit)){
		fit=list()
		rf <- spotPredictRandomForest(rawB,mergedB,largeDesign,spotConfig,NULL)
		mlegp <- spotPredictMlegp(rawB,mergedB,largeDesign,spotConfig,NULL)
		fit[[1]] <-rf$seq.modelFit
		fit[[2]] <-mlegp$seq.modelFit
		newDesignPrediction<-as.data.frame((mlegp$seq.largeDesignY+rf$seq.largeDesignY)/2)
	}else{
		rf.res <- spotPredictRandomForest(NULL,NULL,largeDesign,spotConfig,fit[[1]])$seq.largeDesignY
		mlegp.res <- spotPredictMlegp(NULL,NULL,largeDesign,spotConfig,fit[[2]])$seq.largeDesignY
		newDesignPrediction<-as.data.frame((mlegp.res+rf.res)/2)
	}
	spotConfig$seq.modelFit<-fit;
	spotWriteLines(spotConfig$io.verbosity,3,"spotPredictRandomForestMlegp finished");
	spotConfig$seq.largeDesignY<-newDesignPrediction;
	spotConfig
}