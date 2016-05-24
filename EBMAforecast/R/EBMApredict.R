#' EBMApredict
#'
#' Function allows users to create new predictions given an already estimated EBMA model
#' This function produces predictions based on EBMA model weights and component model predictions.
#'
#' @param EBMAmodel Output of estimated EBMA model 
#' @param Predictions A matrix with a column for each component model's predictions.
#' @param Outcome An optional vector containing the true values of the dependent variable for all observations in the test period. 
#' @param ... Not implemented
#'
#'
#' @return Returns a data of class 'FDatFitLogit' or FDatFitNormal, a subclass of 'ForecastData', with the following slots:
#' \item{predTest}{A matrix containing the predictions of all component models and the EBMA model for all observations in the test period.}#' \item{period}{The period, "calibration" or "test", for which the statistics were calculated.}
#' \item{outcomeTest}{An optional vector containing the true values of the dependent variable for all observations in the test period.}
#' \item{modelNames}{A character vector containing the names of all component models.  If no model names are specified, names will be assigned automatically.}
#' \item{modelWeights}{A vector containing the model weights assigned to each model.}
#' 
#' 
#' @author  Michael D. Ward <\email{michael.d.ward@@duke.edu}> and Jacob M. Montgomery <\email{jacob.montgomery@@wustl.edu}> and Florian M. Hollenbach <\email{florian.hollenbach@@tamu.edu}>
#'
#' @references Montgomery, Jacob M., Florian M. Hollenbach and Michael D. Ward. (2015). Calibrating ensemble forecasting models with sparse data in the social sciences.   \emph{International Journal of Forecasting}. In Press.
#' @references Montgomery, Jacob M., Florian M. Hollenbach and Michael D. Ward. (2012). Improving Predictions Using Ensemble Bayesian Model Averaging. \emph{Political Analysis}. \bold{20}: 271-291.
#'
#' @rdname EBMApredict
#' @aliases prediction, ForecastDataLogit-method prediction, FDatFitNormal-method prediction, FDatFitLogit-method prediction,ForecastDataNormal-method prediction
#' @export
setGeneric(name="EBMApredict",
           def=function(EBMAmodel, 
                        Predictions,
                        Outcome = NULL,
                        ...)
           {standardGeneric("EBMApredict")}
)


#' @rdname EBMApredict
#' @export
setMethod(f="EBMApredict",
          signature="ForecastData",
          definition=function(
            EBMAmodel, 
            Predictions,
            Outcome = NULL,
            ...)
          {
            switch(EBMAmodel@model,
                   logit ={EBMAmodel <- as(EBMAmodel, "FDatFitLogit")},
                   logistic ={EBMAmodel <- as(EBMAmodel, "FDatFitLogit")},
                   normal={EBMAmodel <- as(EBMAmodel, "FDatFitNormal")}
            )
            eval(
              prediction(EBMAmodel,
                      Predictions=Predictions,
                      Outcome=Outcome,
                      ...), parent.frame())
          }
)

