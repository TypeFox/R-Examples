#' ES of Black-Scholes put using Monte Carlo Simulation
#' 
#' Estimates ES of Black-Scholes Put Option using Monte Carlo simulation
#' 
#' @param amountInvested Total amount paid for the Put Option and is positive 
#' (negative) if the option position is long (short)
#' @param stockPrice Stock price of underlying stock
#' @param strike Strike price of the option
#' @param r Risk-free rate
#' @param mu Expected rate of return on the underlying asset and is in 
#' annualised term
#' @param sigma Volatility of the underlying stock and is in annualised 
#' term
#' @param maturity The term to maturity of the option in days
#' @param numberTrials The number of interations in the Monte Carlo simulation
#' exercise
#' @param cl Confidence level for which ES is computed and is scalar
#' @param hp Holding period of the option in days and is scalar
#' @return ES 
#' @references Dowd, Kevin. Measuring Market Risk, Wiley, 2007.
#' 
#' Lyuu, Yuh-Dauh. Financial Engineering & Computation: Principles, 
#' Mathematics, Algorithms, Cambridge University Press, 2002.
#' 
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Market Risk of American Put with given parameters.
#'    BlackScholesPutESSim(0.20, 27.2, 25, .03, .12, .05, 60, 1000, .95, 30)
#'    
#' @export
BlackScholesPutESSim <- function(amountInvested, stockPrice, strike, r, mu, 
                                 sigma, maturity, numberTrials, cl, hp){
  # Precompute Constants
  annualMaturity <- maturity / 360 # Annualised maturity
  annualHp <- hp / 360 # Annualised holding period
  N <- 1 # Number of steps - only one needed for black scholes option
  dt <- annualHp / N # Size of time-increment equal to holding period
  nudt <- (mu - .5 * sigma^2) * dt
  sigmadt <- sigma * sqrt(dt)
  lnS <- log(stockPrice)
  M <- numberTrials
  initialOptionPrice <- BlackScholesPutPrice(stockPrice, strike, 
                                                 r, sigma, maturity)
  numberOfOptions <- abs(amountInvested) / initialOptionPrice
  # Stock price simulation process
  lnSt <- matrix(0, M, N)
  newStockPrice <- matrix(0, M, N)
  for (i in 1:M){
    lnSt[i] <- lnS + rnorm(1, nudt, sigmadt) # Random stock price movement
    newStockPrice[i] <- exp(lnSt[i, 1]) # New stock price
  }
  # Profit/Loss calculation
  profitOrLoss <- double(M)
  if (amountInvested > 0) { # If option position is long
    for (i in 1:M) {
      profitOrLoss[i] <- (BlackScholesPutPrice(newStockPrice[i], strike, r, 
                                              sigma, maturity - hp) - initialOptionPrice) * numberOfOptions
    }
  }
  if (amountInvested < 0) { # If option position is short
    for (i in 1:M) {
      profitOrLoss[i] <- (-BlackScholesPutPrice(newStockPrice[i], strike, r, 
                                              sigma, maturity - hp) + initialOptionPrice) * numberOfOptions
    }
  }
  # VaR estimation
  y <- HSES(profitOrLoss, cl) # VaR
  return(y)
}