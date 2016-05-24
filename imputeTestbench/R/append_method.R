#' To attach and compare new method to existing comparison study done with function 'impute_errors()'
#'
#' @param existing_method as Error observations for different methods
#' @param dataIn as imput time series for testing
#' @param  missPercentFrom as variable from which percent of missing values to be considered
#' @param  missPercentTo as variable to state upto what percent missing values are to be considered
#' @param  interval as interval between consecutive missPercent values
#' @param  repetition as an integer to state number of repetition to be done for each missPercent value
#' @param  errorParameter as type of error calculation (RMSE, MAE or MAPE)
#' @param  MethodPath as location of function for proposed imputation method
#' @param  MethodName as name for function for proposed imputation method
#' @import ggplot2
#' @import imputeTS
#' @importFrom methods hasArg
#' @return Returns error comparosin for imputation methods
#' @export
#' @examples
#' #Kindly, refer "Vignette" document


#==================================================================================
# append_method starts here....
#==================================================================================

append_method <- function(existing_method, dataIn, missPercentFrom, missPercentTo, interval, repetition, errorParameter, MethodPath, MethodName)
{

  if(!(hasArg(dataIn)))
  {
    dataIn <- c(1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5,1:5)
  }

  if(!is.vector(dataIn))
  {
    dataIn <- dataIn[, 1]
  }
  # For future reference
  dataIn1 <- dataIn

  # Set default values
  if(!(hasArg(errorParameter)))
  {
    errorParameter <- 1
  }

  if(!(hasArg(repetition)))
  {
    repetition <- 1
  }


  if(!(hasArg(MethodName)))
  {
    MethodName <- "Proposed Method"
  }

  if(!(hasArg(interval)))
  {
    interval <- 10
  }

  if(!(hasArg(missPercentFrom)))
  {
    missPercentFrom <- 10
  }

  if(!(hasArg(missPercentTo)))
  {
    missPercentTo <- 80
  }

  e <- 0
  f <- 0
  # Function to create missing values
  for(x in seq(missPercentFrom, missPercentTo, interval))
  {
    x <- x/100
    a <- length(dataIn)
    b <- a * x
    b <- abs(b)
    c <- a-b
    out <- NULL

    for(i in 1:repetition)
    {
      dataIn <- dataIn1
      dataIn[c:(c+b)] <- NA
      c <- sample(1:a, 1, replace = TRUE)
      while(c > a-b)
      {
        c <- sample(1:a, 1, replace = TRUE)
      }
      out[i] <- data.frame(dataIn)
    }

    gh <- NULL
    for(i in 1:repetition)
    {
      outs <- as.numeric(unlist(out[i]))

      # to call functions from provided "MethodPath"
        d <- parse(text = MethodPath)
        d <- eval(d)
        d <- d$value(outs)

        if(errorParameter == 1)
        {
          gh[i] <- rmse(dataIn1 - d)
          parameter <- "RMSE Plot"
        }
        if(errorParameter == 2)
        {
          gh[i] <- mae(dataIn1 - d)
          parameter <- "MAE Plot"
        }
        if(errorParameter == 3)
        {
          gh[i] <- mape((dataIn1 - d), dataIn1)
          parameter <- "MAPE Plot"
        }
        if(errorParameter[1] == 4)
        {
          newPar <- parse(text = errorParameter[2])
          newPar <- eval(newPar)
          newPar1 <- newPar$value(dataIn1, d)
          gh[i] <- newPar1
          parameter <- errorParameter[3]
        }
    }



      e <- append(e,mean(gh))
      f <- append(f,x)
  }
    ex <- NULL
    fx <- NULL
    ex <- e[-1]
    fx <- f[-1]
    g <- data.frame(fx,ex)

  #existing_method[length(existing_method)+1] <- data.frame(e[-1])
  existing_method[[paste(MethodName)]] <- e[-1]
  #existing_method[length(existing_method)+1] <- e
  options(warn=-1)
  return(Proposed_Method = existing_method)
}
