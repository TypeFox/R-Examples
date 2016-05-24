#' Calculates the Replacement Cost(RC) and the sum of the MtMs for all the trades
#' @title Calculates the RC
#' @param trades The full list of the Trade Objects 
#' @param coll_agreement (Optional) The collateral Agreement object covering the trade list
#' @param current_collateral (Optional) The current value of the collateral posted from the counterparty to the processing org
#' @return The replacement Cost and the sum of the MtMs
#' @export
#' @author Tasos Grivas <tasos@@openriskcalculator.com>
#' @references Basel Committee: The standardised approach for measuring counterparty credit risk exposures
#' http://www.bis.org/publ/bcbs279.htm

CalcRC <- function(trades,coll_agreement, current_collateral)  {
  
  V <- sum(sapply(trades, function(x) x$MtM))
  
  if (missing(coll_agreement))
  {
    V_C <- V
    RC <- max(V_C,0)
  } else
  {
    V_C <- V - current_collateral
    RC  <- max(V_C,coll_agreement$thres_cpty+coll_agreement$MTA_cpty-coll_agreement$IM_cpty,0)
  }
  
return(list("V_C"=V_C,"RC"=RC))
}