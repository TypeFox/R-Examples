#' Jarque-Bera backtest for normality.
#'
#' Jarque-Bera (JB) is a backtest to test whether the skewness and kurtosis of a
#' given sample matches that of normal distribution. JB test statistic is
#' defined as \deqn{JB=\frac{n}{6}\left(s^2+\frac{(k-3)^2}{4}\right)} where
#' \eqn{n} is sample size, \eqn{s} and \eqn{k} are coefficients of sample
#' skewness and kurtosis.
#' 
#' @param sample.skewness Coefficient of Skewness of the sample
#' @param sample.kurtosis Coefficient of Kurtosis of the sample
#' @param n Number of observations
#' @return Probability of null hypothesis H0
#' 
#' @references Dowd, Kevin. Measuring Market Risk, Wiley, 2007.
#' 
#' Jarque, C. M. and Bera, A. K. A test for normality of observations and
#' regression residuals, International Statistical Review, 55(2): 163-172.
#' 
#' 
#' @author Dinesh Acharya
#' @examples
#' 
#'    # JB test statistic for sample with 500 observations with sample
#'    # skewness and kurtosis of -0.075 and 2.888
#'    JarqueBeraBacktest(-0.075,2.888,500)
#'
#' @export
JarqueBeraBacktest <- function(sample.skewness, sample.kurtosis, n){
  s <- sample.skewness
  k <- sample.kurtosis
  jb.test.stat <- (n/6)*(s^2+((k-3)^2)/4)
  # Corresponding cdf-value for a chi-squared distribution with two degrees of
  # freedom
  prob.value.of.null <- 1-pchisq(jb.test.stat, 2)
  return(prob.value.of.null)
}