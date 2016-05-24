#' VaR for default risky bond portfolio
#'
#' Generates Monte Carlo VaR for default risky bond portfolio in Chapter 6.4
#'  
#' @param r Spot (interest) rate, assumed to be flat
#' @param rf Risk-free rate
#' @param coupon Coupon rate
#' @param sigma Variance
#' @param amount.invested Amount Invested
#' @param recovery.rate Recovery rate
#' @param p Probability of default
#' @param number.trials Number of trials
#' @param hp Holding period
#' @param cl Confidence level
#' @return Monte Carlo VaR
#' 
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#'
#' @author Dinesh Acharya
#' @examples
#' 
#'    # VaR for default risky bond portfolio for given parameters
#'    DefaultRiskyBondVaR(.01, .01, .1, .01, 1, .1, .2, 100, 100, .95)
#'
#' @export
DefaultRiskyBondVaR <- function(r, rf, coupon, sigma, amount.invested, recovery.rate, p, number.trials, hp, cl){
  M <- number.trials
  delta <- recovery.rate
  ann.hp <- hp/360
  # R equals spot rate prevailing at end of hp/2 with term equal to hp.2
  initial.bond.value <- coupon * ((1 + rf)/(1+r))^(ann.hp/2) + (1 + coupon)/((1+r)^ann.hp)
  number.of.bonds <- amount.invested/initial.bond.value
  z <- rnorm(M)
  R <- double(M)
  interim.bond.value <- double(M)
  interim.default.state <- double(M)
  terminal.bond.value <- double(M)
  terminal.default.state <- double(M)
  for (j in 1:M) {
    R[j] <- exp(r + sigma * sqrt(hp/2)*z[j]) # Random realisation of spot rate with term hp/2
    interim.default.state[j] <- rbinom(1,1,p) # Determines whether default occurs at hp/2
    terminal.default.state[j] <- rbinom(1,1,p) # Determines whether default occurs at hp
    if (interim.default.state[j] == 0) {
      if (terminal.default.state[j] == 0) {
        terminal.bond.value[j] <- coupon * (1 + rf) ^ (ann.hp / 2) + (1 + coupon)
      } else {
        terminal.bond.value[j] <- coupon * ((1 + rf) ^ (ann.hp / 2)) + delta * (1 + coupon)
      }
    } else {
      interim.bond.value[j] <- delta * (coupon + ((1 + coupon) / (1+R[j])^(ann.hp/2)))
      terminal.bond.value[j] <- ((1 + rf) ^ (ann.hp / 2)) * interim.bond.value[j]
    }
  }
  profit.or.loss <- number.of.bonds * (terminal.bond.value - initial.bond.value)
  # Convert to P/L
  hist(-profit.or.loss)
  HSVaR(profit.or.loss, cl)
}