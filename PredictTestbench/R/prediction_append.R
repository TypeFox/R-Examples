#' To attach and compare new method to existing comparison study done with function 'impute_errors()'
#'
#' @param existing_method as Error observations for different methods
#' @param dataIn as imput time series for testing
#' @param  nextVal as an integer to decide number of values to predict
#' @param  errorParameter as type of error calculation (RMSE, MAE or MAPE)
#' @param  MethodPath as location of function for proposed imputation method
#' @param  MethodName as name for function for proposed imputation method
#' @import ggplot2
#' @import forecast
#' @import PSF
#' @importFrom imputeTestbench mape
#' @importFrom methods hasArg
#' @return Returns error comparosin for imputation methods
#' @export
#' @examples
#' #Kindly, refer "Vignette" document


#==================================================================================
# prediction_append() starts here....
#==================================================================================

prediction_append <- function(existing_method, dataIn, nextVal, errorParameter, MethodPath, MethodName)
{

  if(!(hasArg(dataIn)))
  {
    dataIn <- c(1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5)
  }

  if(!is.vector(dataIn))
  {
    dataIn <- dataIn[, 1]
  }

  testData <- dataIn[1:(length(dataIn)-nextVal)]
  compData <- dataIn[(length(dataIn)-nextVal+1):length(dataIn)]

  # For future reference
  dataIn1 <- testData

  # Set default values
  if(!(hasArg(errorParameter)))
  {
    errorParameter <- 1
  }

  if(!(hasArg(nextVal)))
  {
    nextVal <- 10
  }

  if(!(hasArg(MethodName)))
  {
    MethodName <- "Proposed Method"
  }


  e <- 0
  f <- 0

    gh <- NULL
    dnewT1 <- Sys.time()
      # to call functions from provided "MethodPath"
      d <- parse(text = MethodPath)
      d <- eval(d)
      d <- d$value(dataIn1)
    dnewT2 <- Sys.time()
    dnewT <- dnewT2 - dnewT1
    dnewT <- as.numeric(dnewT)

      if(errorParameter == 1)
      {
        gh <- rmse(compData - d)
        parameter <- "RMSE Plot"
      }
      if(errorParameter == 2)
      {
        gh <- mae(compData - d)
        parameter <- "MAE Plot"
      }
      if(errorParameter == 3)
      {
        gh <- imputeTestbench::mape((compData - d), compData)
        parameter <- "MAPE Plot"
      }
      if(errorParameter[1] == 4)
      {
        newPar <- parse(text = errorParameter[2])
        newPar <- eval(newPar)
        newPar1 <- newPar$value(compData, d)
        gh <- newPar1
        parameter <- errorParameter[3]
      }

  ex <- NULL
  fx <- NULL
  ex <- d
  fx <- 1:length(compData)
  g <- data.frame(fx,ex)
  MethodName <- paste(MethodName,"Prediction", sep = "_")
  MethodName2 <- paste(MethodName,"Error", sep = "_")
  MethodName3 <- paste(MethodName,"Execution_Time_in_Seconds", sep = "_")

  #existing_method[length(existing_method)+1] <- data.frame(e[-1])
  existing_method[[paste(MethodName)]] <- d
  existing_method[[paste(MethodName2)]] <- gh
  existing_method[[paste(MethodName3)]] <- dnewT

  #existing_method[length(existing_method)+1] <- e
  options(warn=-1)
  return(Proposed_Method = existing_method)
    #return(list(Proposed_Method = existing_method, Proposed_Method_Error = gh, Proposed_Method_Execution_Time_in_Seconds = dnewT))
}
