#' Monte Carlo VaR for DB pension
#' 
#' Generates Monte Carlo VaR for DB pension in Chapter 6.7.
#' 
#' @param mu Expected rate of return on pension-fund assets
#' @param sigma Volatility of rate of return of pension-fund assets
#' @param p Probability of unemployment in any period
#' @param life.expectancy Life expectancy 
#' @param number.trials Number of trials
#' @param cl VaR confidence level
#' @return VaR for DB pension
#' @references Dowd, Kevin. Measuring Market Risk, Wiley, 2007.
#' 
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Estimates the price of an American Put
#'    DBPensionVaR(.06, .2, .05, 80, 100, .95)
#'    
#' @export
DBPensionVaR <- function(mu, sigma, p, life.expectancy, number.trials, cl){
  # Parameter Setting
  contribution.rate <- .15
  initial.income <- 25
  income.growth.rate <- .02
  M <- number.trials
  L <- life.expectancy
  # r is return on investment
  # Asset Side
  # Initialization
  r <- matrix(0, 40, M)
  fund <- matrix(0, 40, M)
  employment.state <- matrix(0, 40, M)
  actual.income <- matrix(0, 40, M)
  contribution <- matrix(0, 40, M)
  years.contributed <- matrix(0, 40, M)
  terminal.fund <- double(M)
  terminal.return <- double(M)
  years.contributed <- matrix(0, 40, M)
  employment.income <-  matrix(0, 40, M)
  total.years.contributed <- double(M)
  for (j in 1:M) {
    fund[1, j] <- contribution.rate * initial.income
    years.contributed[1, j] <- 1
    for (i in 2:40) {
      r[i, j] <- rnorm(1, mu, sigma)
      employment.state[i, j] <- rbinom(1,1,1-p)
      employment.income[i, j] <-  initial.income * exp(income.growth.rate*(i-1))
      actual.income[i, j] <- employment.state[i, j] * employment.income[i, j]
      contribution[i, j] <- contribution.rate * actual.income[i, j]
      fund[i, j] <- contribution[i, j] + fund[i - 1, j] * (1 + r[i, j])
      terminal.fund[j] <- fund[i, j]
      terminal.return[j] <- r[i, j]
      years.contributed[i, j] <- employment.state[i, j] + years.contributed[i - 1, j]
      total.years.contributed[j] <- years.contributed[i,j]
    }
  }
 mean.terminal.fund <- mean(terminal.fund)
 std.terminal.fund <- sd(terminal.fund)
 
 terminal.employment.income <- (1 - p) * initial.income * exp(income.growth.rate * 39)
 pension <- double(M)
 annuity.rate <- double(M)
 implied.fund <- double(M)
 for (j in 1:M){
   pension[j] <- (total.years.contributed[j] / 40) * terminal.employment.income
   annuity.rate[j] <- .04
   implied.fund[j] <- pvfix(annuity.rate[j], L-65, pension[j])
 }
 mean.terminal.employment.income <- mean(terminal.employment.income)
 mean.t.years.contributed <- mean(total.years.contributed)
 mean.pension <- mean(pension)
 mean.implied.func <- mean(implied.fund)
 std.implied.fund <- sd(implied.fund)
 # Profit Loss and VaR
 profit.or.loss <- terminal.fund - implied.fund
 mean.profit.or.loss <- mean(profit.or.loss)
 std.profit.or.loss <- sd(profit.or.loss)
 hist(-profit.or.loss, 20)
 y <- HSVaR(profit.or.loss, cl)
 return(y)
}
# Accessory functions
pvfix<-function(r, n, c){
  # pvfix computes the present value of a series of future cashflows (e.g. savings)
  # parameters:
  # r interest rate per period (constant throughout the period)
  # n number of periods
  # c cashflow each month (assumed to be fixed)
  s <- (c/r)*(1-(1/(1+r)^n))
  return(s)
}