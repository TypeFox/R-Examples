#' VaR of Insurance Portfolio
#'
#' Generates Monte Carlo VaR for insurance portfolio in Chapter 6.5
#'
#' @param mu Mean of returns
#' @param sigma Volatility of returns
#' @param n Number of contracts
#' @param p Probability of any loss event
#' @param theta Expected profit per contract
#' @param deductible Deductible
#' @param number.trials Number of simulation trials
#' @param cl VaR confidence level
#' @return VaR of the specified portfolio
#' 
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#' 
#' 
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Estimates VaR of Insurance portfolio with given parameters
#'    InsuranceVaR(.8, 1.3, 100, .6, 21,  12, 50, .95)
#'
#' @export
InsuranceVaR<- function(mu, sigma, n, p, theta, deductible, number.trials, cl){
  M <- number.trials
  D <- deductible
  L <- matrix(0, n, M)
  company.loss <- L
  total.company.loss <- double(M)
  for (j in 1:M) {
    L[1, j] <- rbinom( 1, 1, p) * rlnorm(1, mu, sigma) # Realisation of L
    company.loss[1, j] <- max(L[1,j] - D, 0) # Adjust for deductible
    
    for (i in 2:n) {
      L[i, j] <- rbinom(1, 1, p) * rlnorm(1, mu, sigma) # Realisation of L
      company.loss[i, j] <- max(L[i,j] - D, 0) + company.loss[i - 1, j] # Adjust
                                                                # for deductible
    }
    total.company.loss[j] <- company.loss[n,j] # Total company loss for 
                                                # given j trial
  }
  # Sample of total company losses
  adjusted.total.company.loss <- total.company.loss - mean(total.company.loss) -
    theta * mean(total.company.loss) / n # Adjusts for premium
  profit.or.loss <- - adjusted.total.company.loss # Convert to P/L
  hist(adjusted.total.company.loss, col = "blue", 
       xlab = "Total Company Loss", ylab = "Frequency",
       main = "Adjusted Total Company Loss")
  y <- HSVaR(profit.or.loss, cl)
  return(y)
  
}