#' Calibrate an ensemble Bayesian Model Averaging model
#'
#' This function calibrates an EBMA model based on out-of-sample performance in the calibration time period. Given a dependent variable and calibration-sample predictions from multiple component forecast models in the \code{ForecastData} the \code{calibrateEnsemble} function fits an ensemble BMA mixture model. The weights assigned to each model are derived from the individual model's performance in the calibration period. Missing observations are allowed in the calibration period, however models with missing observations are penalized. When missing observations are prevalent in the calibration set, the EM algorithm is adjusted and model paprameters are estimated by maximizing a renormalized partial expected complete-data log-likelihood (Fraley et al. 2010).
#'
#' @param .forecastData An object of class 'ForecastData' that will be used to calibrate the model.
#' @param exp The exponential shrinkage term.  Forecasts are raised to the (1/exp) power on the logit scale for the purposes of bias reduction.  The default value is \code{exp=3}.
#' @param tol Tolerance for improvements in the log-likelihood before the EM algorithm will stop optimization.  The default is \code{tol= 0.01}, which is somewhat high.  Researchers may wish to reduce this by an order of magnitude for final model estimation. 
#' @param maxIter The maximum number of iterations the EM algorithm will run before stopping automatically. The default is \code{maxIter=10000}.
#' @param model The model type that should be used given the type of data that is being predicted (i.e., normal, binary, etc.).
#' @param method The estimation method used.  Currently only implements "EM".
#' @param predType The prediction type used for the EBMA model under the normal model, user can choose either \code{posteriorMedian} or \code{posteriorMean}. Posterior median is the default.
#' @param W Vector of initial model weights, if unspecified each model will receive weight 1/number of Models
#' @param const user provided "wisdom of crowds" parameter, serves as minimum model weight for all models. Default = 0
#' @param useModelParams If "TRUE" individual model predictions are transformed based on logit models. If "FALSE" all models' parameters will be set to 0 and 1.  
#' @param ... Not implemented
#'
#' @return Returns a data of class 'FDatFitLogit' or FDatFitNormal, a subclass of 'ForecastData', with the following slots
#' \item{predCalibration}{A matrix containing the predictions of all component models and the EBMA model for all observations in the calibration period.} 
#' \item{predTest}{A matrix containing the predictions of all component models and the EBMA model for all observations in the test period.}
#' \item{outcomeCalibration}{A vector containing the true values of the dependent variable for all observations in the calibration period.} 
#' \item{outcomeTest}{An optional vector containing the true values of the dependent variable for all observations in the test period.}
#' \item{modelNames}{A character vector containing the names of all component models.  If no model names are specified, names will be assigned automatically.}
#' \item{modelWeights}{A vector containing the model weights assigned to each model.}
#' \item{modelParams}{The parameters for the individual logit models that transform the component models.}
#' \item{useModelParams}{Indicator whether model parameters for transformation were estimated or not.}
#' \item{logLik}{The final log-likelihood for the calibrated EBMA model.}
#' \item{exp}{The exponential shrinkage term.}
#' \item{tol}{Tolerance for improvements in the log-likelihood before the EM algorithm will stop optimization.}
#' \item{maxIter}{The maximum number of iterations the EM algorithm will run before stopping automatically.}
#' \item{method}{The estimation method used. }
#' \item{iter}{Number of iterations run in the EM algorithm.}
#' \item{call}{The actual call used to create the object.}
#'
#'
#' @author Michael D. Ward <\email{michael.d.ward@@duke.edu}> and Jacob M. Montgomery <\email{jacob.montgomery@@wustl.edu}> and Florian M. Hollenbach <\email{florian.hollenbach@@tamu.edu}>

#'
#' @references Montgomery, Jacob M., Florian M. Hollenbach and Michael D. Ward. (2015). Calibrating ensemble forecasting models with sparse data in the social sciences.   \emph{International Journal of Forecasting}. In Press.
#' @references Montgomery, Jacob M., Florian M. Hollenbach and Michael D. Ward. (2012). Improving Predictions Using Ensemble Bayesian Model Averaging. \emph{Political Analysis}. \bold{20}: 271-291.
#'
#' @references Raftery, A. E., T. Gneiting, F. Balabdaoui and M. Polakowski. (2005). Using Bayesian Model Averaging to calibrate forecast ensembles. \emph{Monthly Weather Review}. \bold{133}:1155--1174.
#' @references Sloughter, J. M., A. E. Raftery, T. Gneiting and C. Fraley. (2007). Probabilistic quantitative precipitation forecasting using Bayesian model averaging. \emph{Monthly Weather Review}. \bold{135}:3209--3220.
#' @references Fraley, C., A. E. Raftery, T. Gneiting. (2010). Calibrating Multi-Model Forecast Ensembles with Exchangeable and Missing Members using Bayesian Model Averaging. \emph{Monthly Weather Review}. \bold{138}:190--202.
#' @references Sloughter, J. M., T. Gneiting and A. E. Raftery. (2010). Probabilistic wind speed forecasting using ensembles and Bayesian model averaging. \emph{Journal of the American Statistical Association}. \bold{105}:25--35.
#' @references Fraley, C., A. E. Raftery, and T. Gneiting. (2010). Calibrating multimodel forecast ensembles with exchangeable and missing members using Bayesian model averaging. \emph{Monthly Weather Review}. \bold{138}:190--202.
#'
#' @examples \dontrun{data(calibrationSample)
#'
#' data(testSample) 
#'
#' this.ForecastData <- makeForecastData(.predCalibration=calibrationSample[,c("LMER", "SAE", "GLM")],
#'.outcomeCalibration=calibrationSample[,"Insurgency"],.predTest=testSample[,c("LMER", "SAE", "GLM")],
#' .outcomeTest=testSample[,"Insurgency"], .modelNames=c("LMER", "SAE", "GLM"))
#'
#' this.ensemble <- calibrateEnsemble(this.ForecastData, model="logit", tol=0.001, exp=3)
#'}
#'
#' @import abind methods stats graphics
#'
#' @keywords calibrate EBMA 
#'
#' @rdname calibrateEnsemble
#' @aliases fitEnsemble,ForecastData-method fitEnsemble,ForecastDataLogit-method fitEnsemble,ForecastDataNormal-method FDatFitLogit-class ForecastDataLogit-class  ForecastDataNormal-class FDatFitNormal-class calibrateEnsemble,ForecastData-method 
#' @export
setGeneric(name="calibrateEnsemble",
           def=function(.forecastData=new("ForecastData"),
             exp=1,
             tol=sqrt(.Machine$double.eps),
             maxIter=1e6,
             model="logit",
             method="EM",
             ...)
           {standardGeneric("calibrateEnsemble")}
           )



#' @export
setMethod(f="calibrateEnsemble",
          signature="ForecastData",
          definition=function(
            .forecastData,
            exp=1,
            tol=sqrt(.Machine$double.eps),
            maxIter=1e6,
            model="logit",
            method="EM",
            ...)
          {
            switch(model,
                   logit ={.forecastData <- as(.forecastData, "ForecastDataLogit")},
                   logistic ={.forecastData <- as(.forecastData, "ForecastDataLogit")},
                   normal={.forecastData <- as(.forecastData, "ForecastDataNormal")}
                   )
            eval(
              fitEnsemble(.forecastData,
                             exp=exp,
                             tol=tol,
                             maxIter=maxIter,
                             method="EM",
                             ...), parent.frame())
          }
          )




           

