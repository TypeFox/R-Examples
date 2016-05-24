#' Plots cumulative density for KS test and computes confidence interval for
#' KS test stat.
#'
#' Kolmogorov-Smirnov (KS) test statistic is a non parametric test for 
#' distribution equality and measures the maximum distance between two cdfs. 
#' Formally, the KS test statistic is : \deqn{D=\max_i|F(X_i)-\hat{F}(X_i)|}
#' 
#' @param number.trials Number of trials
#' @param sample.size Sizes of the trial samples
#' @param confidence.interval Confidence interval expressed as a fraction of 1
#' @return Confidence Interval for KS test stat
#' 
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#' 
#' Chakravarti, I. M., Laha, R. G. and Roy, J. Handbook of Methods of #' Applied Statistics, Volume 1, Wiley, 1967.
#' 
#' 
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Plots the cdf for KS Test statistic and returns KS confidence interval
#'    # for 100 trials with 1000 sample size and 0.95 confidence interval
#'    KSTestStat(100, 1000, 0.95)
#'
#' @export
KSTestStat <- function(number.trials, sample.size, confidence.interval){
  
  if (confidence.interval>1){
    stop("Confidence Interval should be less than 1.")
  }
  
  # Read back input parameters
  m <- number.trials
  n <- sample.size
  # Random number generator
  data <- matrix(rnorm(m*n), m, n)
  
  # Initialize vectors
  cdf.diff <- double(n)
  max.diff <- double(m)
  
  # Compute KS Test Statistic
  for (i in 1:m) {
    trial.sample <- data[i, ]
    ordered.trial.sample <- sort(trial.sample)
    for (j in 1:n) {
      cdf.diff[j] <- j/n-pnorm(ordered.trial.sample[j],0,1)
    }
    max.diff[i] <- max(abs(cdf.diff))
  }
  max.diff <- sort(max.diff)
  
  # Obtain confidence interval
  lower.bound.index <- round(m*(1-confidence.interval)/2)
  upper.bound.index <- round(m*(confidence.interval+(1-confidence.interval)/2))
  confidence.interval.for.KS.test.stat <- c(max.diff[lower.bound.index], 
                                            max.diff[upper.bound.index])
  
  # Plot
  cdf <- seq(1/m, 1, 1/m)
  plot(max.diff, cdf, col="red", type="l", 
       main="Cumulative Density for KS test statistic", 
       xlab="KS test statistic", ylab="Cumulative probability")
  
  return(confidence.interval.for.KS.test.stat)
}