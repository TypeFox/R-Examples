#' Nonlinearity test
#' @details
#' This function runs a set of nonlinearity tests implemented in other R packages including:
#' \itemize{
#'    \item Teraesvirta's neural metwork test for nonlinearity (\code{\link[tseries]{terasvirta.test}}).
#'    \item White neural metwork test for nonlinearity (\code{\link[tseries]{white.test}}).
#'    \item Keenan's one-degree test for nonlinearity (\code{\link[TSA]{Keenan.test}}).
#'    \item Perform the McLeod-Li test for conditional heteroscedascity (ARCH). (\code{\link[TSA]{McLeod.Li.test}}).
#'    \item Perform the Tsay's test for quadratic nonlinearity in a time series. (\code{\link[TSA]{Tsay.test}}).
#'    \item Perform the Likelihood ratio test for threshold nonlinearity. (\code{\link[TSA]{tlrt}}).
#' }
#' @param time.series The original time.series from which the surrogate data is generated.
#' @param verbose Logical value. If TRUE, a summary of each of the tests is shown.
#' @return A list containing the results of each of the tests.
#' @export nonlinearityTest
#' @import tseries
#' @import TSA
nonlinearityTest <- function(time.series, verbose = TRUE){
  nltests = list()
  # apply all tests
  # testing linearity in mean
  terasvirta = tseries::terasvirta.test(x=ts(time.series),type="Chisq")
  nltests$Terasvirta = terasvirta
  if (verbose){
    cat("\t\t ** Teraesvirta's neural network test  **\n")  
    cat("\t\t Null hypothesis: Linearity in \"mean\" \n")
    cat("\t\t X-squared = ",terasvirta$statistic," df = ",terasvirta$parameter," p-value = ",
        terasvirta$p.value,"\n\n")  
  }  
  white = tseries::white.test(ts(time.series))
  nltests$White = white
  if (verbose){
    cat("\t\t ** White neural network test  **\n")  
    cat("\t\t Null hypothesis: Linearity in \"mean\" \n")
    cat("\t\t X-squared = ",white$statistic," df = ",white$parameter," p-value = ",
        white$p.value,"\n\n")    
  }
  # nonlinearity against the null hypothesis that the time series follows some AR process.
  keenan = TSA::Keenan.test(time.series)
  nltests$Keenan = keenan
  if (verbose){
    cat("\t\t ** Keenan's one-degree test for nonlinearity  **\n")  
    cat("\t\t Null hypothesis: The time series follows some AR process\n")
    cat("\t\t F-stat = ",keenan$test.stat," p-value = ", keenan$p.value,"\n\n")    
  }
  # test for conditional heteroscedascity
  mcleodLi = TSA::McLeod.Li.test(y=time.series,plot=FALSE) 
  nltests$McLeodLi = mcleodLi
  if (verbose){
    cat("\t\t ** McLeod-Li test  **\n")  
    cat("\t\t Null hypothesis: The time series follows some ARIMA process\n")
    cat("\t\t Maximum p-value = ",  max( unlist(mcleodLi) ),"\n\n")    
  }
  # Tsay test for nonlinearity
  tsay = TSA::Tsay.test(time.series)
  nltests$Tsay = tsay
  if (verbose){
    cat("\t\t ** Tsay's Test for nonlinearity **\n")  
    cat("\t\t Null hypothesis: The time series follows some AR process\n")
    cat("\t\t F-stat = ",tsay$test.stat," p-value = ", tsay$p.value,"\n\n")    
  }
  # Likelihood ratio test for threshold nonlinearity
  tarTest = TSA::tlrt(time.series)
  nltests$TarTest = tarTest
  if (verbose){
    cat("\t\t ** Likelihood ratio test for threshold nonlinearity **\n")  
    cat("\t\t Null hypothesis: The time series follows some AR process\n")
    cat("\t\t Alternativce hypothesis: The time series follows some TAR process\n")
    cat("\t\t X-squared = ",tarTest$test.stat," p-value = ", tarTest$p.value,"\n\n")    
  }
  
  
  return(nltests)
  
}

