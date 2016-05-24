#' @description Loads the supervisory data (factors, correlation and option volatility)
#' for each Asset Class and SubClass
#' @title Supervisory Data Loading
#'  
#' @return A data frame with the required data
#' @export
#' @author Tasos Grivas <tasos@@openriskcalculator.com>
#' @references Basel Committee: The standardised approach for measuring counterparty credit risk exposures
#' http://www.bis.org/publ/bcbs279.htm
LoadSupervisoryData <- function()  {
  
   # reading the excel file containing the supervisory factors
  superv <- read.csv(system.file("extdata", "supervisory_factors.csv", package = "SACCR"),header=TRUE,stringsAsFactors = FALSE)
  superv$Supervisory_factor = as.numeric(sub("%","",superv$Supervisory_factor))/100
  superv$Correlation = as.numeric(sub("%","",superv$Correlation))/100
  superv$Supervisory_option_volatility = as.numeric(sub("%","",superv$Supervisory_option_volatility))/100
  return(superv)
  
}