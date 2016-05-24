#' Binomial Put Price
#' 
#' Estimates the price of an American Put, using the binomial approach.
#' 
#' @param stockPrice Stock price of underlying stock
#' @param strike Strike price of the option
#' @param r Risk-free rate
#' @param sigma Volatility of the underlying stock and is in annualised 
#' term
#' @param maturity The term to maturity of the option in days
#' @param numberSteps The number of time-steps in the binomial tree
#' @return Binomial American put price
#' @references Dowd, Kevin. Measuring Market Risk, Wiley, 2007.
#' 
#' Lyuu, Yuh-Dauh. Financial Engineering & Computation: Principles, 
#' Mathematics, Algorithms, Cambridge University Press, 2002.
#' 
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Estimates the price of an American Put
#'    AmericanPutPriceBinomial(27.2, 25, .03, .2, 60, 30)
#'    
#' @export
AmericanPutPriceBinomial <- function(stockPrice, strike, r, sigma, 
                             maturity, numberSteps){
  prob <- .5 # Probability of up-move, equal probability approach used
  N <- numberSteps
  maturity <- maturity/360 # Convert maturity to units of years
  deltat <- maturity/N # Length of incremental steps
  u <- 2*exp(r * deltat + 2 * sigma * sqrt(deltat)) / (exp(2 * sigma * sqrt(deltat)) + 1) # Up-move
  d <- 2*exp(r * deltat) / (exp(2 * sigma * sqrt(deltat)) + 1) # Dowd-move
  # Note that there up- and down- moves incorporate risk  neutral approach
  discount <- exp(-r * deltat) # Discount factor
  # Stock price tree
  S <- double(N+1)
  S[1] <- stockPrice * (d ^ N)
  for (i in 2:(N + 1)) { # i is number of up-moves plus 1
    S[i] <- S[i-1] * u * u
  }
  # Option price tree
  # Terminal option values
  x <- double(N + 1)
  for (i in 1:(N + 1)) {
    x[i] <- max(strike - S[i], 0)
  }
  # Option values prior to expiration
  for (j in seq(N, 1, -1)) {
    for (i in 1:j){
      x[i] <- discount * (prob * x[i + 1] + (1 - prob) * x[i])
      S[i] <- S[i] * u
      # Check for early exercise
      exerciseValue <- max(0, strike - S[i])
      if (exerciseValue > x[i]) {
        x[i] <- exerciseValue
      }
    }
  }
  y <- x[1] # Initial option value
  return(y)
}