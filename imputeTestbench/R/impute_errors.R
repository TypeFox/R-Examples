#' Function working as testbench for comparison of imputing models
#'
#' @param  dataIn as imput time series for testing
#' @param  missPercentFrom as variable from which percent of missing values to be considered
#' @param  missPercentTo as variable to state upto what percent missing values are to be considered
#' @param  interval as interval between consecutive missPercent values
#' @param  repetition as an integer to decide the numbers of repetition to be done for each missPercent value
#' @param  errorParameter as type of error calculation (RMSE, MAE or MAPE)
#' @param  MethodPath as location of function for the proposed imputation method
#' @param  MethodName as name for function for the proposed imputation method
#' @import ggplot2
#' @import imputeTS
#' @importFrom stats ts
#' @importFrom methods hasArg
#' @return Returns error comparison for imputation methods
#' @export
#' @examples
#' # aa <- impute_errors()
#' # aa

#==================================================================================
# impute_error starts here....
#==================================================================================

impute_errors <- function(dataIn, missPercentFrom, missPercentTo, interval, repetition, errorParameter, MethodPath, MethodName)
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
  e1 <- 0
  f1 <- 0
  enew <- 0
  fnew <- 0

  # Function to create missing values
  for(x in seq(missPercentFrom, missPercentTo, interval))
  {
    x <- x/100
    # Inputs: dataIn, miss_per, repetition
    #repetition <- 5      # number of repetition
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
    gh1 <- NULL
    ghnew <- NULL
    for(i in 1:repetition)
    {
      outs <- as.numeric(unlist(out[i]))

      #d <- impute(outs,mean)
      out <- ts(outs)
      d <- na.mean(out)
      d1 <- na.interpolation(out)

      if((hasArg(MethodPath)))
      {
        # to call functions from provided "MethodPath"
        dnew <- parse(text = MethodPath)
        dnew <- eval(dnew)
        dnew <- dnew$value(outs)

        if(errorParameter == 1)
        {
          ghnew[i] <- rmse(dataIn1 - dnew)
          parameter <- "RMSE Plot"
        }
        if(errorParameter == 2)
        {
          ghnew[i] <- mae(dataIn1 - dnew)
          parameter <- "MAE Plot"
        }
        if(errorParameter == 3)
        {
          ghnew[i] <- mape((dataIn1 - dnew), dataIn1)
          parameter <- "MAPE Plot"
        }
        if(errorParameter[1] == 4)
        {
          newPar <- parse(text = errorParameter[2])
          newPar <- eval(newPar)
          newPar <- newPar$value(dataIn1, dnew)
          ghnew[i] <- newPar
          parameter <- errorParameter[3]
        }
      }


      if(errorParameter == 1)
      {
        gh[i] <- rmse(dataIn1 - d)
        gh1[i] <- rmse(dataIn1 - d1)
        parameter <- "RMSE Plot"
      }
      if(errorParameter == 2)
      {
        gh[i] <- mae(dataIn1 - d)
        gh1[i] <- mae(dataIn1 - d1)
        parameter <- "MAE Plot"
      }
      if(errorParameter == 3)
      {
        gh[i] <- mape((dataIn1 - d), dataIn1)
        gh1[i] <- mape((dataIn1 - d1), dataIn1)
        parameter <- "MAPE Plot"
      }
      if(errorParameter[1] == 4)
      {
        newPar <- parse(text = errorParameter[2])
        newPar <- eval(newPar)
        newPar1 <- newPar$value(dataIn1, d)
        gh[i] <- newPar1
        newPar2 <- newPar$value(dataIn1, d1)
        gh1[i] <- newPar2
        parameter <- errorParameter[3]
      }
    }

    #gh
    e <- append(e,mean(gh))
    f <- append(f,x)
    e1 <- append(e1,mean(gh1))
    f1 <- append(f1,x)

    #x <- x + 10
    if((hasArg(MethodPath)))
    {
      enew <- append(enew,mean(ghnew))
      fnew <- append(fnew,x)
    }
  }
  e
  ex <- e[-1]
  fx <- f[-1]
  g <- data.frame(fx,ex)
  #plot(e, type = 'l')
  h <- ggplot(g,aes(fx,ex)) + labs(title = parameter) + xlab("Percent of Missing Values")+ ylab("Error Values") + geom_line(aes(color="Historic Mean")) + labs(color="Imputing Methods")

  ex <- e1[-1]
  fx <- f1[-1]
  g1 <- data.frame(fx,ex)
  h <- h + geom_line(data=g1,aes(color= "Interpolation"))

  if((hasArg(MethodPath)))
  {
    ex <- NULL
    fx <- NULL
    ex <- enew[-1]
    fx <- fnew[-1]
    gnew <- data.frame(fx,ex)

    h <- h + geom_line(data=gnew,aes(color= MethodName))
  }
  options(warn=-1)

  #axis(1, labels = e)
 #return(list(Missing_Percent = fx[-1], Historic_Mean = e[-1], Interpolation = e1[-1], Proposed_Method = enew[-1], Plot = h))
  if((hasArg(MethodPath)))
  {

  return(list(Parameter = parameter, Missing_Percent = f[-1], Historic_Mean = e[-1], Interpolation = e1[-1], Proposed_Method = enew[-1]))
  }
  else{
    return(list(Parameter = parameter, Missing_Percent = f[-1], Historic_Mean = e[-1], Interpolation = e1[-1]))
  }
 }
