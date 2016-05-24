#' Calculates the xVA values (CVA, DVA, FVA, FBA, MVA, KVA)
#'
#' @title Calculates the xVA values
#' @param trades The full list of the Trade Objects
#' @param col    The margin agreement with the counterparty
#' @param sim_data A list containing data related to the calculation of simulated exposures (for example the model parameters and the number of simulations)
#' @param reg_data A list containing data related to the regulatory calculations (for example the 'framework' member variable can be 'IMM','SACCR','CEM')
#' @param credit_curve_PO   The credit curve of the processing organisation
#' @param credit_curve_cpty The credit curve of the processing organisation
#' @param funding_curve     A curve containing the credit spread for the funding of the collateral
#' @param spot_rates        The spot rates curve
#' @param cpty_LGD          The loss-given-default of the counterparty
#' @param PO_LGD            The loss-given-default of the processing organisation
#' @return A list containing the xVA values
#' @export
#' @author Tasos Grivas <tasos@@openriskcalculator.com>
#' @references Gregory J., The xVA Challenge, 2015, Wiley
#'
xVACalculator = function(trades, col, sim_data, reg_data, credit_curve_PO, credit_curve_cpty, funding_curve, spot_rates, cpty_LGD, PO_LGD)
{
  maturity      <- max(as.numeric(lapply(trades, function(x) x$Ei)))
  time_points    = GenerateTimeGrid(col, maturity)
  num_of_points  = length(time_points)

  spot_curve     = spot_rates$CalcInterpPoints(time_points)
  cpty_spread    = credit_curve_cpty$CalcInterpPoints(time_points)
  PO_spread      = credit_curve_PO$CalcInterpPoints(time_points)
  funding_spread = funding_curve$CalcInterpPoints(time_points)

  discount_factors = exp(-time_points*spot_curve)

  PD_cpty        = CalcPD(cpty_spread,cpty_LGD,time_points)
  PD_PO          = CalcPD(PO_spread,PO_LGD,time_points)
  PD_FVA         = CalcPD(funding_spread,1,time_points)

  exposure_profile = CalcSimulatedExposure(discount_factors, time_points, spot_curve, col, trades, sim_data)

  xVA = list()

  xVA$KVA        = calcKVA(exposure_profile, col,trades, reg_data, time_points)
  xVA$CVA        = CalcVA(exposure_profile$EE,  discount_factors, PD_cpty, cpty_LGD)
  xVA$DVA        = CalcVA(exposure_profile$NEE, discount_factors, PD_PO, PO_LGD)
  xVA$FCA        = CalcVA(exposure_profile$EE,  discount_factors, PD_FVA)
  xVA$FBA        = CalcVA(exposure_profile$NEE, discount_factors, PD_FVA)
  xVA$MVA        = xVA$FCA*2*sqrt(reg_data$mva_days/(250*maturity))*qnorm(reg_data$mva_percentile)/dnorm(0)

  return(xVA)
}
