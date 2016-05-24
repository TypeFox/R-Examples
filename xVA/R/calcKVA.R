#' Calculates the capital valuation adjustment by computing the default capital charge and the CVA capital charge and applying the required return-on-capital
#' @title Calculates the Capital Valuation Adjustment (KVA)
#' @param exposure_profile The exposure profile list containing the EE, EEE etc
#' @param col    The margin agreement with the counterparty
#' @param trades The full list of the Trade Objects
#' @param reg_data A list containing data related to the regulatory calculations (for example the 'framework' member variable can be 'IMM','SACCR','CEM')
#' @param time_points The timepoints that the analysis is performed on
#' @return The capital valuation adjustment (KVA)
#' @export
#' @author Tasos Grivas <tasos@@openriskcalculator.com>
#'
calcKVA = function(exposure_profile, col,trades, reg_data, time_points)
{

  EAD = calcEAD(trades, reg_data$framework, col,exposure_profile$EEE,time_points)

  effective_maturity = calcEffectiveMaturity(trades, time_points, reg_data$framework, exposure_profile$EE)

  def_capital_charge = calcDefCapital(trades,EAD, reg_data, effective_maturity)

  cva_capital_charge = calcCVACapital(trades, EAD, reg_data$cpty_rating, effective_maturity)

  KVA = -(def_capital_charge+cva_capital_charge)*0.5*sqrt(effective_maturity)*reg_data$return_on_capital
}
