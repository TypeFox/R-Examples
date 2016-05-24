#' Calculates the probablity of the default on specific time points by using the spread of the corresponding credit curve and the loss given default
#' @title Calculates the Probablity of Default
#' @param spread The spread based on the credit curve
#' @param LGD The loss-given-default
#' @param time_points The timepoints that the analysis is performed on
#' @return A vector containing the probablity of default on the specified timepoints
#' @export
#' @author Tasos Grivas <tasos@@openriskcalculator.com>
#'

CalcPD = function(spread, LGD, time_points)
{
  num_of_points = length(time_points)
  spread = spread/10^4
  PD = exp(-(spread[1:(num_of_points-1)]*time_points[1:(num_of_points-1)])/LGD) - exp(-(spread[2:num_of_points]*time_points[2:num_of_points])/LGD)

  PD[PD<0] = 0
  return(PD)
}
