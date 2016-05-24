#' Calculates the Exposure at Default
#' @title Calculates the EAD
#' @param RC        the replacement cost
#' @param PFE       the projected future exposure
#' @return The Exposure-at-Default
#' @export
#' @examples
#' #returns 1.4*(60+500) = 784
#' EAD <- CalcEAD(60,500)
#' @author Tasos Grivas <tasos@@openriskcalculator.com>
#' @references Basel Committee: The standardised approach for measuring counterparty credit risk exposures
#' http://www.bis.org/publ/bcbs279.htm
CalcEAD <- function(RC, PFE)  {
  
  EAD <- 1.4*(RC+PFE)
  
  return(EAD)
}
