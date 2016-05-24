#' Function working as testbench for comparison of Prediction algorithms
#'
#' @param  dataIn as input time series for testing
#' @param  nextVal as an integer to decide number of values to predict
#' @param  errorParameter as type of error calculation (RMSE, MAE or MAPE)
#' @param  MethodPath as location of function for the proposed imputation method
#' @param  MethodName as name for function for the proposed imputation method
#' @import ggplot2
#' @import forecast
#' @import PSF
#' @importFrom imputeTestbench mape
#' @importFrom stats ts
#' @importFrom methods hasArg
#' @return Returns error comparison for imputation methods
#' @export
#' @examples
#' # aa <- prediction_errors(nextVal = 10)
#' # aa

#==================================================================================
# prediction_errors starts here....
#==================================================================================

prediction_errors <- function(dataIn, nextVal, errorParameter, MethodPath, MethodName)
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
  e1 <- 0
  f1 <- 0
  enew <- 0
  fnew <- 0

    gh <- NULL
    gh1 <- NULL
    ghnew <- NULL

    dT1 <- Sys.time()
      d <- forecast(auto.arima(dataIn1))
      d <- as.numeric(unlist(data.frame(d)[1]))
    dT2 <- Sys.time()
    dT <- dT2 - dT1
    dT <- as.numeric(dT)

    d1T1 <- Sys.time()
      d1 <- AUTO_PSF(dataIn1, nextVal)$Predicted_Values
    d1T2 <- Sys.time()
    d1T <- d1T2 - d1T1
    d1T <- as.numeric(d1T)

     if((hasArg(MethodPath)))
      {
        # to call functions from provided "MethodPath"
    dnewT1 <- Sys.time()
        dnew <- parse(text = MethodPath)
        dnew <- eval(dnew)
        dnew <- dnew$value(dataIn1)
    dnewT2 <- Sys.time()
    dnewT <- dnewT2 - dnewT1
    dnewT <- as.numeric(dnewT)

        if(errorParameter == 1)
        {
          ghnew <- rmse(compData - dnew)
          parameter <- "RMSE Plot"
        }
        if(errorParameter == 2)
        {
          ghnew <- mae(compData - dnew)
          parameter <- "MAE Plot"
        }
        if(errorParameter == 3)
        {
          ghnew <- imputeTestbench::mape((compData - dnew), compData)
          parameter <- "MAPE Plot"
        }
        if(errorParameter[1] == 4)
        {
          newPar <- parse(text = errorParameter[2])
          newPar <- eval(newPar)
          newPar <- newPar$value(compData, dnew)
          ghnew <- newPar
          parameter <- errorParameter[3]
        }
      }


      if(errorParameter == 1)
      {
        gh <- rmse(compData - d)
        gh1 <- rmse(compData - d1)
        parameter <- "RMSE Plot"
      }
      if(errorParameter == 2)
      {
        gh <- mae(compData - d)
        gh1 <- mae(compData - d1)
        parameter <- "MAE Plot"
      }
      if(errorParameter == 3)
      {
        gh <- imputeTestbench::mape((compData - d), compData)
        gh1 <- imputeTestbench::mape((compData - d1), compData)
        parameter <- "MAPE Plot"
      }
      if(errorParameter[1] == 4)
      {
        newPar <- parse(text = errorParameter[2])
        newPar <- eval(newPar)
        newPar1 <- newPar$value(compData, d)
        gh <- newPar1
        newPar2 <- newPar$value(compData, d1)
        gh1 <- newPar2
        parameter <- errorParameter[3]
      }
  e
  fx <- 1:length(compData)
  ex <- d
  g <- data.frame(fx,ex)
  h <- ggplot(g,aes(fx,ex)) + labs(title = parameter) + xlab("Percent of Missing Values")+ ylab("Error Values") + geom_line(aes(color="ARIMA Method")) + labs(color="Imputing Methods")

  ex <- d1
  g1 <- data.frame(fx,ex)
  h <- h + geom_line(data=g1,aes(color= "PSF Method"))

  if((hasArg(MethodPath)))
  {
    ex <- enew
    fx <- 1:length(compData)
    gnew <- data.frame(fx,ex)

    h <- h + geom_line(data=gnew,aes(color= MethodName))
  }
  options(warn=-1)

  if((hasArg(MethodPath)))
  {

    return(list(Parameter = parameter, Desired_Prediction = compData, ARIMA_Method_Prediction = d, ARIMA_Method_Error = gh, ARIMA_Execution_Time_in_Seconds = dT, PSF_Method_Prediction = d1, PSF_Method_Error = gh1, PSF_Execution_Time_in_Seconds = d1T, Proposed_Method_Prediction = dnew, Proposed_Method_Error = ghnew, Proposed_Method_Execution_Time_in_Seconds = dnewT))
  }
  else{
    return(list(Parameter = parameter, Desired_Prediction = compData, ARIMA_Method_Prediction = d, ARIMA_Method_Error = gh, ARIMA_Execution_Time_in_Seconds = dT, PSF_Method_Prediction = d1, PSF_Method_Error = gh1, PSF_Execution_Time_in_Seconds = d1T))
    }
}
