#' Calculates the CVA capital charge based on the standardized approach
#' @title Calculates the CVA Capital Charge
#' @param trades The full list of the Trade Objects
#' @param EAD Exposure-at-Default
#' @param cpty_rating the rating of the counterparty
#' @param effective_maturity The effective maturity of the trades of the netting set
#' @return The CVA capital charge of the trade set
#' @export
#' @author Tasos Grivas <tasos@@openriskcalculator.com>
calcCVACapital = function(trades, EAD, cpty_rating, effective_maturity)
{
  rating_table = HashTable('RatingsMapping.csv',"character","numeric")

  reg_weight =rating_table$FindValue(cpty_rating)

  df =(1-exp(-0.05*effective_maturity))/(effective_maturity*0.05)

  cva_capital_charge = qnorm(0.99)*sqrt((0.5*reg_weight*EAD*df*effective_maturity)^2+0.75*(reg_weight*EAD*df*effective_maturity)^2)

  return(cva_capital_charge)
}
