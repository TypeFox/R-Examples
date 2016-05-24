#' Calculates the Valuation Adjustment based on the exposure, the probability-of-default and the loss-given-default
#' @title Calculates the Valuation Adjustment
#' @param exposure A vector containing the exposure values on which the credit risk adjustment will be calculated
#' @param discount_factors The Discount Curve
#' @param PD       The probability-of-Default
#' @param LGD      The Loss-Given-Default
#' @return The Valuation Adjustment Value
#' @export
#' @author Tasos Grivas <tasos@@openriskcalculator.com>
#'
CalcVA = function(exposure, discount_factors, PD, LGD)
{
  VA = 0
  num_of_points = length(discount_factors)

  if(missing(LGD))
    LGD = 1

  for (i in 1:(num_of_points-1))
    VA = VA + exposure[i]*discount_factors[i]*PD[i]

  VA = -VA*LGD

}
