#' Build a ensemble forecasting data object
#'
#' This function uses the component model forecasts and dependent variable observations provided by the user to create an object of class \code{ForecastData}, which can then be used to calibrate and fit the ensemble. Individual slots of the \code{ForecastData} object can be accessed and changed using the \code{get} and \code{set} functions respectively. Missing predictions are allowed in the calibration set.
#'
#' @param .predCalibration A matrix with the number of rows being the number of observations in the calibration period and a column with calibration period predictions for each model.
#' @param .predTest A vector with the number of rows being the number of observations in the test period and a column with test period predictions for each model.
#' @param .outcomeCalibration A vector with the true values of the dependent variable for each observation in the calibration period.  
#' @param .outcomeTest A vector with the true values of the dependent variable for each observation in the test period.
#' @param .modelNames A vector of length p with the names of the component models.  
#' @param ... Additional arguments not implemented
#'
#' Additionally, the functions \code{show} and \code{print} can be used to display data objects of class 'ForecastData'.
#' \code{show} displays only 1 digit and takes the following parameters:
#' @param x A data object of class 'ForecastData'
#' 
#' \code{print} let's the use specify the number of digits printed and takes the arguments:
#' @param object A data object of class 'ForecastData'
#' @param digits User specified number of digits to be displayed.
#'
#' @return A data object of the class 'ForecastData' with the following slots: 
#' \item{predCalibration}{An array containing the predictions of all component models for all observations in the calibration period.} 
#' \item{predTest}{An array containing the predictions of all component models for all observations in the test period.}
#' \item{outcomeCalibration}{A vector containing the true values of the dependent variable for all observations in the calibration period.} 
#' \item{outcomeTest}{A vector containing the true values of the dependent variable for all observations in the test period.}
#' \item{modelNames}{A character vector containing the names of all component models.  If no model names are specified, names will be assigned automatically.}
#' @author  Michael D. Ward <\email{michael.d.ward@@duke.edu}> and Jacob M. Montgomery <\email{jacob.montgomery@@wustl.edu}> and Florian M. Hollenbach <\email{florian.hollenbach@@tamu.edu}>
#'
#'
#' @references Montgomery, Jacob M., Florian M. Hollenbach and Michael D. Ward. (2015). Calibrating ensemble forecasting models with sparse data in the social sciences.   \emph{International Journal of Forecasting}. In Press.
#' @references Montgomery, Jacob M., Florian M. Hollenbach and Michael D. Ward. (2012). Improving Predictions Using Ensemble Bayesian Model Averaging. \emph{Political Analysis}. \bold{20}: 271-291.
#'
#' @examples  data(calibrationSample)
#' 
#' \dontrun{data(testSample) 
#' this.ForecastData <- makeForecastData(.predCalibration=calibrationSample[,c("LMER", "SAE", "GLM")],
#' .outcomeCalibration=calibrationSample[,"Insurgency"],.predTest=testSample[,c("LMER", "SAE", "GLM")],
#' .outcomeTest=testSample[,"Insurgency"], .modelNames=c("LMER", "SAE", "GLM"))
#' 
#' ### to acces individual slots in the ForecastData object
#' getPredCalibration(this.ForecastData)
#' getOutcomeCalibration(this.ForecastData)
#' getPredTest(this.ForecastData)
#' getOutcomeTest(this.ForecastData)
#' getModelNames(this.ForecastData)
#' 
#' ### to assign individual slots, use set functions
#'
#' setPredCalibration(this.ForecastData)<-calibrationSample[,c("LMER", "SAE", "GLM")]
#' setOutcomeCalibration(this.ForecastData)<-calibrationSample[,"Insurgency"]
#' setPredTest(this.ForecastData)<-testSample[,c("LMER", "SAE", "GLM")]
#' setOutcomeTest(this.ForecastData)<-testSample[,"Insurgency"]
#' setModelNames(this.ForecastData)<-c("LMER", "SAE", "GLM")
#'}
#'
#' @seealso ensembleBMA
#' @aliases makeForecastData-method, ForecastData-method setModelNames<-, ForecastData-generic setModelNames<-, setOutcomeTest<-, setOutcomeCalibration, setPredTest, setPredCalibration, print-method, show-method
#' @rdname makeForecastData
#' @export
setGeneric(name="makeForecastData",
           def=function(
            .predCalibration=array(NA, dim=c(0,0,0)),
             .predTest=array(NA, dim=c(0,0,0)),
             .outcomeCalibration=numeric(),
            .outcomeTest=numeric(),
             .modelNames=character(),
             ...)
           {standardGeneric("makeForecastData")}
           )


#' @rdname makeForecastData
#' @export
setMethod(f="makeForecastData",
          definition=function(
            .predCalibration,
            .predTest,
            .outcomeCalibration,
            .outcomeTest,
            .modelNames)
          {
            if(class(.predCalibration)=="data.frame"){.predCalibration <- as.matrix(.predCalibration)}
            if(class(.predTest)=="data.frame"){.predTest <- as.matrix(.predTest)}
            if(class(.predCalibration)=="matrix"){.predCalibration <- array(.predCalibration, dim=c(nrow(.predCalibration), ncol(.predCalibration), 1))}
            if(class(.predTest)=="matrix"){.predTest <- array(.predTest, dim=c(nrow(.predTest), ncol(.predTest), 1))}
            if(length(.modelNames)<ncol(.predCalibration)){
              .modelNames <- paste("Model", 1:ncol(.predCalibration))
            }
            if(length(.predCalibration)>0){
              colnames(.predCalibration) <- .modelNames; rownames(.predCalibration) <- 1:nrow(.predCalibration)
            }
            if (length(.predTest)>0){
              colnames(.predTest) <- .modelNames; rownames(.predTest) <- 1:nrow(.predTest)
            }
            if(length(.outcomeCalibration>0)) {names(.outcomeCalibration) <- 1:length(.outcomeCalibration)}
            if(length(.outcomeTest>0))  {names(.outcomeTest) <- 1:length(.outcomeTest)}
            
            return(new("ForecastData", predCalibration=.predCalibration, predTest=.predTest,
                       outcomeCalibration=.outcomeCalibration, outcomeTest=.outcomeTest, modelNames=.modelNames))
            
          }
          )

#' @rdname makeForecastData
#' @export
setMethod(
		f="print",
		signature="ForecastData",
		definition=function(x, digits=3, ...){
			cat("* Prediction Calibration = \n"); 
			if(length(x@predCalibration)>0)
			{print(x@predCalibration, na.print="", digits=digits);}
			else{print("Nothing Here")}
			cat("* Prediction Test = \n"); 
			if(length(x@predTest)>0)
			{print(x@predTest, na.print="", digits=digits);}
			else{print("Nothing Here")}
			cat("* Outcome Calibration = \n");
			if(length(x@outcomeCalibration)>0)
			{print(x@outcomeCalibration, na.print="", digits=digits);}
			else{print("Nothing Here")}
			cat("* Outcome Test = \n");
			if(length(x@outcomeTest)>0)
			{print(x@outcomeTest, na.print="", digits=digits);}
			else{print("Nothing Here")}
			cat("* Model Names = \n ");print(x@modelNames, na.print="");
			}
			)

#' @rdname makeForecastData
#' @export
setMethod(
		f="show",
		signature="ForecastData",
		definition=function(object){
                  if (length(object@predCalibration)==0) {
			cat("* Prediction Calibration = \n"); 
			if(length(object@predCalibration)>0)
			{print(object@predCalibration, na.print="", digits=1);}
			else{print("Nothing Here")}
			cat("* Prediction Test = \n"); 
			if(length(object@predTest)>0)
			{print(object@predTest, na.print="", digits=1);}
			else{print("Nothing Here")}
			cat("* Outcome Calibration = \n");
			if(length(object@outcomeCalibration)>0)
			{print(object@outcomeCalibration, na.print="", digits=1);}
			else{print("Nothing Here")}
			cat("* Outcome Test = \n");
			if(length(object@outcomeTest)>0)
			{print(object@outcomeTest, na.print="", digits=1);}
			else{print("Nothing Here")}
			cat("* Model Names = \n ");print(object@modelNames, na.print="");
                  }
            else{
            nrowCal=min(10,nrow(object@predCalibration))
            nrowTest=min(10,nrow(object@predTest))
            	cat("* Prediction Calibration = \n"); 
				if(length(object@predCalibration)>0)
				{print(object@predCalibration[1:nrowCal,1:ncol(object@predCalibration),1], na.print="", digits=2);}
				else{print("Nothing Here")}
				cat("* Prediction Test = \n"); 
				if(length(object@predTest)>0)
				{print(object@predTest[1:nrowTest,1:ncol(object@predTest),1], na.print="", digits=2);}
				else{print("Nothing Here")}
				cat("* Outcome Calibration = \n");
				if(length(object@outcomeCalibration)>0)
				{print(print(object@outcomeCalibration[1:nrowCal]),na.print="", digits=2);}
				else{print("Nothing Here")}
				cat("* Outcome Test = \n");
				if(length(object@outcomeTest)>0)
				{print(object@outcomeTest[1:nrowTest], na.print="", digits=2);}
				else{print("Nothing Here")}
				cat("* Model Names = \n ");print(object@modelNames,na.print="");
            	}
            }
          )
