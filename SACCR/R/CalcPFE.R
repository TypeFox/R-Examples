#' Calculates the Projected Future Exposure (PFE) after applying the relevant multiplier.
#' The purpose of this multiplier is to lessen the risk stemming from the addons in case of excess collateral
#' @title Calculates the PFE
#' @param V_C the difference between the sum of the MtMs and the collateral
#' @param Addon_Aggregate the aggregate amount of the Addon
#' @return The Projected Future Exposure (PFE)
#' @export
#' @author Project team <info@@openriskcalculator.com>
#' @references Basel Committee: The standardised approach for measuring counterparty credit risk exposures
#' http://www.bis.org/publ/bcbs279.htm

CalcPFE <- function(V_C,Addon_Aggregate)  {
  
  if (V_C<0){
    multiplier <- min( 1, 0.05 + 0.95 * exp(V_C/(1.9*Addon_Aggregate)))
  }else  
  { multiplier <- 1}
  
  PFE <- multiplier * Addon_Aggregate
  
  return(PFE)
}