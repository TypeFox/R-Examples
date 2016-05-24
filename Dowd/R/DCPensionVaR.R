#' Monte Carlo VaR for DC pension
#' 
#' Generates Monte Carlo VaR for DC pension in Chapter 6.7.
#' 
#' @param mu Expected rate of return on pension-fund assets
#' @param sigma Volatility of rate of return of pension-fund assets
#' @param p Probability of unemployment in any period
#' @param life.expectancy Life expectancy 
#' @param number.trials Number of trials
#' @param cl VaR confidence level
#' @return VaR for DC pension
#' @references Dowd, Kevin. Measuring Market Risk, Wiley, 2007.
#' 
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Estimates the price of an American Put
#'    DCPensionVaR(.06, .2, .05, 80, 100, .95)
#'    
#' @export
DCPensionVaR <- function(mu, sigma, p, life.expectancy, number.trials, cl){
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
  annuity.rate <- double(M)
  pension <- double(M)
  pension.ratio <- double(M)
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
      annuity.rate[j] <- .04
      pension[j] <- payper(annuity.rate[j], L - 65, terminal.fund[j])
      terminal.employment.income <- (1 - p) * initial.income * exp(income.growth.rate * 39)
      pension.ratio[j] <- pension[j]/terminal.employment.income
    }
  }
  mean.terminal.fund <- mean(terminal.fund)
  std.terminal.fund <- sd(terminal.fund)
  # Histogram
  hist(pension.ratio, 20)
  # VaR
  VaR <- -HSVaR(pension.ratio, cl)
  return(VaR)
}

# Access Function
payper <- function(r, n, p){
  # Computes payment per period for annuity or loans
  # parameters:
  # r interest rate per period
  # n number of periods
  # p present value of the instrument
  payment <- (p * r) / (1-(1/(1+r)^n))
  return(payment)
}