#' Calculates the effective maturity based on the specified regulatory framework
#' @title Calculates the Effective Maturity
#' @param trades The full list of the Trade Objects
#' @param framework Specifies the regulatory framework used in the calculations. It can take the values of 'IMM', 'CEM', 'SA-CCR'
#' @param simulated_exposure The exposure profile list containing the EE, EEE etc
#' @param time_points The timepoints that the analysis is performed on
#' @return The effective maturity of the trade set
#' @export
#' @author Tasos Grivas <tasos@@openriskcalculator.com>
#'
calcEffectiveMaturity = function(trades, time_points, framework, simulated_exposure)
{
  if(framework == "IMM")
  {
    effective_maturity = max(1,min( sum(simulated_exposure[time_points>1]*diff(time_points[time_points>=1]))/sum(simulated_exposure[time_points<1]*diff(time_points[time_points<=1])),5))

  }
  else
  {
    Notional_vector = as.numeric(lapply(trades, function(x) x$Notional))

    if(framework == "SA-CCR")
    {
      adj_notional_vector = as.numeric(lapply(trades, function(x) x$CalcAdjNotional()))
      effective_maturity = max(1, min(sum(adj_notional_vector)/sum(Notional_vector),5))
    }
    else if(framework == "CEM")
    {
      Maturity_vector = as.numeric(lapply(trades, function(x) x$Ei))
      effective_maturity = max(1, min(5,sum(Notional_vector*Maturity_vector)/sum(Notional_vector)))
    }
  }

  return(effective_maturity)
}
