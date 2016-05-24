#' Log Normal VaR with stop loss limit
#'
#' Generates Monte Carlo lognormal VaR with stop-loss limit
#'  
#' @param mu Mean arithmetic return
#' @param sigma Standard deviation of arithmetic return
#' @param number.trials Number of trials used in the simulations
#' @param loss.limit Stop Loss limit
#' @param cl Confidence Level
#' @param hp Holding Period
#' @return Lognormal VaR
#' 
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#'
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Estimates standard error of normal quantile estimate
#'    StopLossLogNormalVaR(0, .2, 100, 1.2, .95, 10)
#'
#' @export
StopLossLogNormalVaR <- function(mu, sigma, number.trials, loss.limit, cl, hp){
  N <- 100 # Number of increments, taken as 100
  dt <- hp/N # Size of time-increment
  nudt <- (mu - 0.5 * sigma ^ 2) * dt
  sigmadt <- sigma * sqrt(dt)
  stock.price <- 1 # Stock price assumed to be investment assumed to be 1
  lnS <- log(stock.price)
  M <- number.trials
  L <- loss.limit
  # Stock price simulation process
  lnSt <- matrix(0,M,N)
  new.stock.price <- matrix(0,M,N)
  equity.proportion <- matrix(0,M,N)
  investment <- matrix(0,M,N)
  for (i in 1:M) {
    lnSt[i, 1] <- rnorm(1, lnS + nudt, sigmadt)
    new.stock.price[i,1] <- exp(lnSt[i,1])
    for (j in 2:N) {
      lnSt[i, j] <- rnorm(1, lnSt[i, j-1] + nudt, sigmadt)
      new.stock.price[i, j] <- exp(lnSt[i, j]) # New stock price
      investment[i, j] <- max(new.stock.price[i, j], stock.price - L)
    }
  }
  # Profit/Loss calculation
  profit.or.loss <- double(M)
  for (i in 1:M) {
    profit.or.loss[i] <- investment[i,j] - stock.price
  }
  y <- HSVaR(profit.or.loss, cl)
  return(y)
}