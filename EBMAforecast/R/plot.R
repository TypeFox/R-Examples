##
#' Plotting function for ensemble models of the class "FDatFitLogit" or "FDatFitNormal", which are the objects created by the \code{calibrateEnsemble()} function.
#'
#' Default plotting for objectes created by the "calibrateEnsemble" function.  See details below.
#'
#' For objects of the class "FDatFitLogit", this function creates separation plots for each of the fitted models, including the EBMA model. Observations are ordered from left to right with increasing predicted probabilities, which is depicted by the black line. Actual occurrences are displayed by red vertical lines. Plots can be displayed for the test or calibration period.
#' For objects of the class "FDatFitNormal", this function creates a plot of the predictive density distribution containing the EBMA PDF and the PDFs for all component models (scaled by their model weights).  It also plots the prediction for the ensemble and the components for the specified observations.
#'
#' @param x An object of class "FDatFitLogit" or "FDatFitNormal"
#' @param period Can take value of "calibration" or "test" and indicates the period for which the plots should be produced.
#' @param subset The row names or numbers for the observations the user wishes to plot.  Only implemented for the subclass "FDatFitNormal"
#' @param mainLabel A vector strings to appear at the top of each predictive posterior plot.  Only implemented for the subclass "FDatFitNormal"
#' @param xLab The label for the x-axis. Only implemented for the subclass "FDatFitNormal"
#' @param yLab The label for the y-axis.  Only implemented for the subclass "FDatFitNormal"
#' @param cols A vector containing the color for plotting the predictive pdf of each component model forecast. Only implemented for the subclass "FDatFitNormal" 
#' @param ... Not implemented
#' 
#' @method plot FDatFitLogit
#' @method plot FDatFitNormal
#' @return NULL
#'
#' @author  Michael D. Ward <\email{michael.d.ward@@duke.edu}> and Jacob M. Montgomery <\email{jacob.montgomery@@wustl.edu}> and Florian M. Hollenbach <\email{florian.hollenbach@@tamu.edu}>
#'
#' @references Raftery, A. E., T. Gneiting, F. Balabdaoui and M. Polakowski. (2005). Using Bayesian Model Averaging to calibrate forecast ensembles. \emph{Monthly Weather Review}. \bold{133}:1155--1174.
#' @references Greenhill, B., M.D. Ward, A. Sacks. (2011). The Separation Plot: A New Visual Method For Evaluating the Fit of Binary Data. \emph{American Journal of Political Science}.\bold{55}: 991--1002.
#' @references Montgomery, Jacob M., Florian M. Hollenbach and Michael D. Ward. (2012). Improving Predictions Using Ensemble Bayesian Model Averaging. \emph{Political Analysis}. \bold{20}: 271-291.
#' @references Montgomery, Jacob M., Florian M. Hollenbach and Michael D. Ward. (2015). Calibrating ensemble forecasting models with sparse data in the social sciences.   \emph{International Journal of Forecasting}. In Press.
#'
#' @seealso \code{separationplot}
#'
#' @import separationplot
#'
#' @examples data(calibrationSample)
#'
#' data(testSample) 
#' 
#' this.ForecastData <- makeForecastData(.predCalibration=calibrationSample[,c("LMER", "SAE", "GLM")],
#' .outcomeCalibration=calibrationSample[,"Insurgency"],.predTest=testSample[,c("LMER", "SAE", "GLM")],
#' .outcomeTest=testSample[,"Insurgency"], .modelNames=c("LMER", "SAE", "GLM"))
#' 
#' this.ensemble <- calibrateEnsemble(this.ForecastData, model="logit", tol=0.001, exp=3)
#' 
#' plot(this.ensemble, period="calibration") 
#' plot(this.ensemble, period="test")
#'
#' @aliases plot,FDatFitLogit-method plot,FDatFitNormal-method
#' @export
setMethod(
          f="plot",
          signature="FDatFitLogit",
          definition=function(x, period="calibration",  subset=1, 
          mainLabel="", xLab="", yLab="", cols=1, ...){ #everything behind period, is not used for logit and is a hack so we don't get a warning about unused documented functions
            nDraw=1
            numModels <- length(x@modelWeights)+1
            modelNames <- c("EBMA", x@modelNames)
            if(period=="calibration"){
              .pred <- x@predCalibration; .actual <- x@outcomeCalibration
            }
            else{
              .pred <- x@predTest; .actual <- x@outcomeTest
            }
            par(mgp=c(1, 0, 0), lend = 2, mar=c(1,0,1,0), mfrow=c(numModels, 1))
            for (i in 1:numModels){
              .miss <- is.na(.pred[,i, nDraw])
              separationplot(pred=as.vector(.pred[!.miss,i, nDraw]), actual=as.vector(.actual[!.miss]), heading=modelNames[i], newplot=F)
            }
          }
          )




#' @export
setMethod(
          f="plot",
          signature="FDatFitNormal",
          definition=function(x, period="calibration",  subset=1,
            mainLabel=paste("Observation", subset), xLab="Outcome", yLab="Posterior Probability", cols=2:(length(x@modelNames)+1), ... )
          {

            thisDraw=1

            if(period=="calibration"){
              .nMod <- length(x@modelNames)
              .pred <- matrix(x@predCalibration[subset,,thisDraw], ncol=.nMod+1);  colnames(.pred) <- c("EBMA", x@modelNames)
              .actual <- x@outcomeCalibration[subset]
            } else{
              .nMod <- length(x@modelNames)
              .pred <- matrix(x@predTest[subset,,thisDraw], ncol=.nMod+1); colnames(.pred) <- c("EBMA", x@modelNames)
              .actual <- x@outcomeTest
            }
            
            .sd <- sqrt(x@variance)
            
            if (length(subset)>1){
              for (j in 1:nrow(.pred)){
                .means <- .pred[j,x@modelNames]
                .miss <- is.na(.means)
                .nModThis <- sum(!.miss)
                .means <- .means[!.miss]
                
                .xMin <- min(.means)-2.5*.sd;  .xMax <- max(.means)+2.5*.sd
                .xAxis <- seq(.xMin, .xMax, by=.01);  .yAxis <- matrix(NA, .nModThis, length(.xAxis)) 
                W <- x@modelWeights[!.miss]
                for(i in 1:.nModThis){ .yAxis[i,] <- dnorm(.xAxis, mean=.means[i], sd=.sd)*W[i] }
                .totals <- colSums(.yAxis)
                plot(NULL, xlim=c(.xMin, .xMax), ylim=c(0,max(.totals)), main=mainLabel[j], xlab=xLab, ylab=yLab)
                for(i in 1:.nModThis){lines(.xAxis, .yAxis[i,], type="l", lty=2,  col=cols[i])}
                lines(.xAxis, colSums(.yAxis), lwd=2)
                rug(.means);  rug(.pred[j,"EBMA"], lwd=3)
                abline(v=.actual[j], lty=3)
              }
            } else {
              .means <- .pred[,x@modelNames]
              .miss <- is.na(.means)
              .nModThis <- sum(!.miss)
              .means <- .means[!.miss]
              .xMin <- min(.means)-2.5*.sd
              .xMax <- max(.means)+2.5*.sd
              .xAxis <- seq(.xMin, .xMax, by=.01)
              .yAxis <- matrix(NA, .nModThis, length(.xAxis)) 
              W <- x@modelWeights[!.miss]
              for(i in 1:.nModThis){.yAxis[i,] <- dnorm(.xAxis, mean=.means[i], sd=.sd)*W[i]}
              .totals <- colSums(.yAxis)
              plot(NULL, xlim=c(.xMin, .xMax), ylim=c(0,max(.totals)), main=mainLabel, xlab=xLab, ylab=yLab)
              for(i in 1:.nModThis){lines(.xAxis, .yAxis[i,], type="l", lty=2, col=cols[i])}
              lines(.xAxis, colSums(.yAxis))
              rug(.means); rug(.pred[,"EBMA"], lwd=3)
              abline(v=.actual, lty=3)
            }
          }
          )



