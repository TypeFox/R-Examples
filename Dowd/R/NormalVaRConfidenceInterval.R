#' Generates Monte Carlo 95\% Confidence Intervals for normal VaR
#' 
#' Generates 95\% confidence intervals for normal VaR using Monte Carlo simulation
#' 
#' @param mu Mean of the P/L process
#' @param sigma Standard deviation of the P/L process
#' @param number.trials Number of trials used in the simulations
#' @param sample.size Sample drawn in each trial
#' @param cl Confidence Level
#' @param hp Holding Period
#'
#' @return 95\% confidence intervals for normal VaR
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#' 
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Generates 95\% confidence intervals for normal VaR for given parameters
#'    NormalVaRConfidenceInterval(0, .5, 20, 15, .95, 90)
#'    
#' 
#' @export
NormalVaRConfidenceInterval <- function(mu, sigma, number.trials, sample.size, cl, hp){
  VaR <- double(number.trials)
  for (k in 1:number.trials) {
    z <- rnorm(sample.size)
    x <- sigma * z + mu
    VaR[k] <- NormalVaR(returns = x, cl = cl, hp  = hp)
  }
  VaR <- sort(VaR)
  lower.order.stat <- floor(0.025 * number.trials)
  upper.order.stat <- ceiling(0.975 * number.trials)
  y <- c(VaR[lower.order.stat], VaR[upper.order.stat])
  return(y)
}
