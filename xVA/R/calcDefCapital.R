#' Calculates the default capital charge using the advanced IRB methodology and the stressed R
#' @title Calculates the Default Capital Charge
#' @param trades The full list of the Trade Objects
#' @param reg_data A list containing data related to the regulatory calculations (for example the regulatory probability-of-default, the regulatory loss-given-default etc)
#' @param effective_maturity          The effective maturity of the trades of the netting set
#' @param EAD                         The Exposure-At-Default of the trades as per the selected regulatory framework
#' @return The default capital charge
#' @export
#' @author Tasos Grivas <tasos@@openriskcalculator.com>
#'
calcDefCapital = function(trades,EAD, reg_data, effective_maturity)
{
  R_bII  = 0.12*(1-exp(-50*reg_data$PD))/(1-exp(-50))+0.24*(1-(1-exp(-50*reg_data$PD))/(1-exp(-50)))
  R_stressed = 1.25*R_bII
  b = (0.11852-0.05478*log(reg_data$PD))^2
  K = (pnorm((qnorm(reg_data$PD)/sqrt(1-R_stressed) + sqrt(R_stressed/(1-R_stressed))*qnorm(0.999)))-reg_data$PD)*(1+(effective_maturity-2.5)*b)/(1-1.5*b)

  RWA = 12.5*K*EAD

  default_capital_charge = 0.08*reg_data$LGD*RWA

  return(default_capital_charge)
}
