#' An ensemble forecasting data object
#'
#' Objects of class \code{ForecastData} are used in the \code{calibrateEnsemble} function. Datasets should be converted into an object of class \code{ForecastData} using the \code{makeForecastData} function. Individual slots of the \code{ForecastData} object can be accessed and changed using the \code{get} and \code{set} functions respectively. Missing observations in the prediction calibration set are allowed.
#'
#'
#' A data object of the class 'ForecastData' has the following slots: 
#'  @slot predCalibration An array containing the predictions of all component models for the observations in the calibration period.
#'  @slot predTest An array containing the predictions of all component models for all the observations in the test period.
#'  @slot outcomeCalibration A vector containing the true values of the dependent variable for all observations in the calibration period. 
#'  @slot outcomeTest A vector containing the true values of the dependent variable for all observations in the test period.
#'  @slot modelNames A character vector containing the names of all component models. 
#'
#' @author  Michael D. Ward <\email{michael.d.ward@@duke.edu}> and Jacob M. Montgomery <\email{jacob.montgomery@@wustl.edu}> and Florian M. Hollenbach <\email{florian.hollenbach@@tamu.edu}>  
#'
#' @references Montgomery, Jacob M., Florian M. Hollenbach and Michael D. Ward. (2015). Calibrating ensemble forecasting models with sparse data in the social sciences.   \emph{International Journal of Forecasting}. In Press.
#' @references Montgomery, Jacob M., Florian M. Hollenbach and Michael D. Ward. (2012). Improving Predictions Using Ensemble Bayesian Model Averaging. \emph{Political Analysis}. \bold{20}: 271-291.
#' 
#' @examples \dontrun{ data(calibrationSample)
#' 
#' data(testSample) 
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
#' @param object used for validity checks (internal)
#' @param value used for validity checks (internal) #hack but no idea how to get the warning to go away otherwise
#'
#'
#' @docType class
#' @seealso ensembleBMA
#' @rdname ForecastData  




setClass(Class="ForecastData",
         representation = representation(
           predCalibration="array",
           predTest="array",
           outcomeCalibration="numeric",
           outcomeTest="numeric",
           modelNames="character"),
         prototype=prototype(
           predCalibration=array(NA, dim=c(0,0,0)),
           predTest=array(NA, dim=c(0,0,0)), 
           outcomeCalibration=numeric(),
           outcomeTest=numeric(),
           modelNames=character()),
         validity=function(object){
         	if(length(object@predCalibration)>0 | length(object@predTest)>0 ){
           		if(nrow(object@predCalibration)!=length(object@outcomeCalibration))
             	{stop("The number of predictions and outcomes do not match in the calibration set.")}
           	}
            #if(length(object@predTest)>0 | length(object@outcomeTest)>0){ 
           	#	if(nrow(object@predTest)!=length(object@outcomeTest))
            # 	{warning("The number of predictions and outcomes do not match in the test set.", call.=FALSE)}
           	#} 
           	if(length(object@predTest)>0 & length(object@predCalibration)>0){
             	if(ncol(object@predTest)!=ncol(object@predCalibration))
               	{stop("The number of prediction models in the calibration and test set are different.")}    
             	if(dim(object@predTest)[3]!=dim(object@predCalibration)[3])
               	{stop("The number of exchangeable draws per model in the calibration and test are different.")}
           	}
           	if(sum(is.na(object@outcomeCalibration)) > 0){
             {stop("There are NAs in the outcome calibration set, these observations should be deleted from the data.")}
           	if(sum(is.na(object@outcomeTest)) > 0)
             {stop("There are NAs in the outcome test set, these observations should be deleted from the data.")}
            }
            if(any(apply(object@predCalibration,c(1,3),FUN=function(x){all(is.na(x)==TRUE)}))==TRUE){
              {stop("One of the observations in the calibration set has missing values for all prediction models.")}
            }
            if(any(apply(object@predTest,c(1,3),function(x){all(is.na(x)==TRUE)}))==TRUE){
              {stop("One of the observations in the test set has missing values for all prediction models.")}
            }
         }
)


##
setMethod("initialize", "ForecastData", function(.Object, ...) {
  value = callNextMethod()
  validObject(value)
  return(value)
})


setClass(Class="ForecastDataLogit",
         contains="ForecastData",
         validity=function(object){
           if(any(object@outcomeCalibration!=1 & object@outcomeCalibration!=0 & !is.na(object@outcomeCalibration))) 
             {stop("The outcomes for the binary model should be either 0 or 1 (Not true for outcome calibration set).")}	
           if(any(object@outcomeTest!=1 & object@outcomeTest!=0 & !is.na(object@outcomeTest))) 
             {stop("The outcomes for the binary model should be either 0 or 1 (Not true for outcome test set).")}	
           if(any(object@predCalibration<0 & !is.na(object@predCalibration)) | (any(object@predCalibration>1 & !is.na(object@predCalibration))))                 {stop("The predictions for the binary model should be between 0 or 1 (Not true for prediction calibration set).")}	
           if(any(object@predTest<0 & !is.na(object@predTest)) |any(object@predTest>1 & !is.na(object@predTest)) )
             {stop("The predictions for the binary model should be between 0 or 1 (Not true for prediction test set).")}	
           if(any(object@predCalibration==0, na.rm=TRUE) | any(object@predCalibration==1, na.rm=TRUE)) 
             {stop("The predictions for the binary model cannot be exactly 0 or 1 (Not true for prediction calibration set).")}	
           if(any(object@predTest==0, na.rm=TRUE) | any(object@predTest==1, na.rm=TRUE)) 
             {stop("The predictions for the binary model cannot be exactly 0 or 1 (Not true for prediction test set).")}	
           
                 }
         )

setClass(Class="ForecastDataNormal",
         contains="ForecastData")


##
setAs(from="ForecastData", to="ForecastDataLogit",
      def=function(from){
        new("ForecastDataLogit",
            predCalibration=from@predCalibration,
            predTest=from@predTest,
            outcomeCalibration=from@outcomeCalibration,
            outcomeTest=from@outcomeTest,
            modelNames=from@modelNames)
}
)

##
setAs(from="ForecastData", to="ForecastDataNormal",
      def=function(from){
        new("ForecastDataNormal",
            predCalibration=from@predCalibration,
            predTest=from@predTest,
            outcomeCalibration=from@outcomeCalibration,
            outcomeTest=from@outcomeTest,
            modelNames=from@modelNames)
      }
      )

