#' Estimates ES of American vanilla put using binomial option valuation tree and Monte Carlo
#' Simulation
#' 
#' Estimates ES of American Put Option using binomial tree to price the option
#' valuation tree and Monte Carlo simulation with a binomial option valuation 
#' tree nested within the MCS. Historical method to compute the VaR.
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
#' @param numberSteps The number of steps over the holding period at each
#' of which early exercise is checked and is at least 2
#' @param cl Confidence level for which VaR is computed and is scalar
#' @param hp Holding period of the option in days and is scalar
#' @return Monte Carlo Simulation VaR estimate and the bounds of the 95% 
#' confidence interval for the VaR, based on an order-statistics analysis 
#' of the P/L distribution
#' @references Dowd, Kevin. Measuring Market Risk, Wiley, 2007.
#' 
#' Lyuu, Yuh-Dauh. Financial Engineering & Computation: Principles, 
#' Mathematics, Algorithms, Cambridge University Press, 2002.
#' 
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Market Risk of American Put with given parameters.
#'    AmericanPutESSim(0.20, 27.2, 25, .16, .2, .05, 60, 30, 20, .95, 30)
#'    
#' @export
AmericanPutESSim <- function(amountInvested, stockPrice, strike, r, mu, sigma, 
                              maturity, numberTrials, numberSteps, cl, hp){
  # Precompute Constants
  annualMaturity <- maturity / 360 # Annualised maturity
  annualHp <- hp / 360 # Annualised holding period
  N <- numberSteps # Number of steps
  dt <- annualHp / N # Size of time-increment equal to holding period
  nudt <- (mu - .5 * sigma^2) * dt
  sigmadt <- sigma * sqrt(dt)
  lnS <- log(stockPrice)
  M <- numberTrials
  initialOptionPrice <- AmericanPutPriceBinomial(stockPrice, strike, r, sigma, maturity, N)
  numberOfOptions <- abs(amountInvested) / initialOptionPrice
  # Stock price simulation process
  lnSt <- matrix(0, M, N)
  newStockPrice <- matrix(0, M, N)
  for (i in 1:M){
    lnSt[i, 1] <- lnS + rnorm(1, nudt, sigmadt)
    newStockPrice[i, 1] <- exp(lnSt[i, 1])
    for (j in 2:N){
      lnSt[i, j] <- lnSt[i, j - 1] + rnorm(1, nudt, sigmadt)
      newStockPrice[i, j] <- exp(lnSt[i,j]) # New stock price
    }
  }
  
  # Option calculation over time
  newOptionValue <- matrix(0, M, N-1)
  for (i in 1:M) {
    for (j in 1:(N-1)) {
      newOptionValue[i, j] <- AmericanPutPriceBinomial(newStockPrice[i, j], 
                                                       strike, r, sigma, maturity - j * hp / N, N)
    }
  }
  # Profit/Loss
  profitOrLoss <- (newOptionValue - initialOptionPrice)*numberOfOptions
  
  # Now adjust for short position
  if (amountInvested < 0) {# If option position is short
    profitOrLoss <- -profitOrLoss
  }
  
  # VaR estimation
  ES <- HSESDFPerc(profitOrLoss, .5, cl) # VaR
  confidenceInterval <- c(HSESDFPerc(profitOrLoss, .025, cl), HSESDFPerc(profitOrLoss, .975, cl))
  return(list('ES' = ES, 'confidenceInterval' = confidenceInterval))
}