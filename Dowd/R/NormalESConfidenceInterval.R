#' Generates Monte Carlo 95\% Confidence Intervals for normal ES
#' 
#' Generates 95\% confidence intervals for normal ES using Monte Carlo simulation
#' 
#' @param mu Mean of the P/L process
#' @param sigma Standard deviation of the P/L process
#' @param number.trials Number of trials used in the simulations
#' @param sample.size Sample drawn in each trial
#' @param cl Confidence Level
#' @param hp Holding Period
#'
#' @return 95\% confidence intervals for normal ES
#' 
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#' 
#' @author Dinesh Acharya
#' 
#' @examples
#' 
#'    # Generates 95\% confidence intervals for normal ES for given parameters
#'    NormalESConfidenceInterval(0, .5, 20, 20, .95, 90)
#'    
#' 
#' @export
NormalESConfidenceInterval <- function(mu, sigma, number.trials, sample.size, cl, hp){
  ES <- double(number.trials)
  for (k in 1:number.trials) {
    z <- rnorm(sample.size)
    x <- sigma * z + mu
    ES[k] <- NormalES(returns = x, cl = cl, hp = hp)
  }
  ES <- sort(ES)
  lower.order.stat <- floor(0.025 * number.trials)
  upper.order.stat <- ceiling(0.975 * number.trials)
  y <- c(ES[lower.order.stat], ES[upper.order.stat])
  return(y)
}
