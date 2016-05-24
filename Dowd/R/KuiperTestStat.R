#' Plots cummulative density for Kuiper test and computes confidence interval 
#' for Kuiper test stat.
#'
#' Kuiper test statistic is a non parametric test for 
#' distribution equality and is closely related to KS test. Formally, the 
#' Kuiper test statistic is :
#'  \deqn{D*=\max_i\{F(X_i)-\hat{F(x_i)}+\max_i\{\hat{F}(X_i)-F(X_i)\}}
#' 
#' @param number.trials Number of trials
#' @param sample.size Sizes of the trial samples
#' @param confidence.interval Confidence interval expressed as a fraction of 1
#' @return Confidence Interval for KS test stat
#' 
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#' 
#' 
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Plots the cdf for Kuiper Test statistic and returns Kuiper confidence 
#'    # interval for 100 trials with 1000 sample size and 0.95 confidence 
#'    # interval.
#'    KuiperTestStat(100, 1000, 0.95)
#'
#' @export
KuiperTestStat <- function(number.trials, sample.size, confidence.interval){
  
  if (confidence.interval >= 1) {
    stop("Confidence Interval should be less than 1.")
  }
  if (confidence.interval <= 0) {
    stop("Confidence Interval should be positive.")
  }
  
  # Read back input parameters
  m <- number.trials
  n <- sample.size
  
  # Random number generator
  data <- matrix(rnorm(m*n), m, n)
  
  # Initialize vectors
  cdf.diff <- double(n)
  Kuiper.test.stat <- double(m)
  
  # Compute KS Test Statistic
  for (i in 1:m) {
    trial.sample <- data[i, ]
    ordered.trial.sample <- sort(trial.sample)
    for (j in 1:n) {
      cdf.diff[j] <- j/n-pnorm(ordered.trial.sample[j],0,1)
    }
    Kuiper.test.stat[i] <- max(cdf.diff)-min(cdf.diff)
  }
  Kuiper.test.stat <- sort(Kuiper.test.stat)
  
  # Obtain confidence interval
  lower.bound.index <- round(m*(1-confidence.interval)/2)
  upper.bound.index <- round(m*(confidence.interval+(1-confidence.interval)/2))
  confidence.interval.for.Kuiper.test.stat <- c(Kuiper.test.stat[lower.bound.index], 
                                            Kuiper.test.stat[upper.bound.index])
  
  # Plot
  cdf <- seq(1/m, 1, 1/m)
  plot(Kuiper.test.stat, cdf, col="red", type="l", 
       main="Cumulative density for Kuiper test statistic", 
       xlab="Kuiper test statistic", ylab="Cumulative probability")
  
  return(confidence.interval.for.Kuiper.test.stat)
}