#' Plots cumulative density for AD test and computes confidence 
#' interval for AD test stat.
#' 
#' Anderson-Darling(AD) test can be used to carry out distribution equality test and is 
#' similar to Kolmogorov-Smirnov test. AD test statistic is defined as:
#' \deqn{A^2=n\int_{-\infty}^{\infty}\frac{[\hat{F}(x)-F(x)]^2}{F(x)[1-F(x)]}dF(x)}
#' which is equivalent to
#' \deqn{=-n-\frac{1}{n}\sum_{i=1}^n(2i-1)[\ln F(X_i)+\ln(1-F(X_{n+1-i}))]}
#' 
#' @param number.trials Number of trials
#' @param sample.size Sample size
#' @param confidence.interval Confidence Interval
#' @return Confidence Interval for AD test statistic
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#' 
#' Anderson, T.W. and Darling, D.A. Asymptotic Theory of Certain Goodness of
#' Fit Criteria Based on Stochastic Processes, The Annals of Mathematical
#' Statistics, 23(2), 1952, p. 193-212.
#' 
#' Kvam, P.H. and Vidakovic, B. Nonparametric Statistics with Applications to
#' Science and Engineering, Wiley, 2007.
#' 
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Probability that the VaR model is correct for 3 failures, 100 number
#'    # observations and  95% confidence level
#'    ADTestStat(1000, 100, 0.95)
#'
#' @export
ADTestStat <- function(number.trials, sample.size, confidence.interval){
  
  if (confidence.interval >= 1){
    stop("Confidence Interval should be less than 1.")
  }
  if (confidence.interval <= 0){
    stop("Confidence Interval should be positive.")
  }
  
  m <- number.trials
  n <- sample.size
  
  # Random Number Generation
  data <- matrix(rnorm(m*n), m, n)
  
  # Initialize vectors
  term <- double(n)
  AD.test.stat <- double(m)
  
  # Compute AD test statistic
  for (i in 1:m) {
    trial.sample <- data[i, ]
    ordered.trial.sample <- sort(trial.sample)
    for (j in 1:n) {
      term[j] <- (2*j-1)*(log(pnorm(ordered.trial.sample[j],0,1))-
                            log(1-pnorm(ordered.trial.sample[n+1-j],0,1)));
    }
    AD.test.stat[i] <- -n-mean(term)
  }
  AD.test.stat <- sort(AD.test.stat)
  
  # Obtain confidence interval
  lower.bound.index <- round(m*(1-confidence.interval)/2)
  upper.bound.index <- round(m* (confidence.interval+(1-confidence.interval)/2))
  confidence.interval.for.AD.test.stat <- c(AD.test.stat[lower.bound.index], 
                                            AD.test.stat[upper.bound.index])
  # Plot the graph
  cdf <- seq(1/m, 1, 1/m)
  plot(AD.test.stat, cdf, col="red", type="l", 
       main="Cumulative density for AD test statistic", 
       xlab="AD test statistic", ylab="Cumulative probability")
  
  return(confidence.interval.for.AD.test.stat)
  
}
