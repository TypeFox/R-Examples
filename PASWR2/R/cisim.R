#' @title Confidence Interval Simulation Program
#' 
#' @description This program simulates random samples from which it constructs confidence intervals for either the population mean, the population variance, or the population proportion of successes. 
#' 
#' @details Default is to construct confidence intervals for the population mean. Simulated confidence intervals for the population variance or population proportion of successes are possible by selecting the appropriate value in the \code{type} argument.
#' 
#' @param samples the number of samples desired.
#' @param n the size of each sample
#' @param parameter If constructing confidence intervals for the population mean or the population variance, parameter is the population mean (i.e., type is one of either \code{"Mean"} or \code{"Var"}). If constructing confidence intervals for the population proportion of successes, the value entered for parameter represents the population proportion of successes (\code{Pi}), and as such, must be a number between 0 and 1.
#' @param sigma is the population standard deviation. \code{sigma} is not required if confidence intervals are of type \code{"Pi"}.
#' @param conf.level confidence level for the graphed confidence intervals, restricted to lie between zero and one
#' @param type character string, one of \code{"Mean"}, \code{"Var"}, or \code{"Pi"}, or just the initial letter of each, indicating the type of confidence interval simulation to perform
#' 
#' @return Performs specified simulation and draws the resulting confidence intervals on a graphical device.
#' 
#' @author Alan T. Arnholt <arnholtat@@appstate.edu> 
#' 
#' @export
#' 
#' @examples
#' cisim(samples = 100, n = 30, parameter = 100, sigma = 10, conf.level = 0.90)
#' # Simulates 100 samples of size 30 from  a normal distribution with mean 100
#' # and a standard deviation of 10.  From the 100 simulated samples, 90% confidence
#' # intervals for the Mean are constructed and depicted in the graph. 
#' 
#' cisim(100, 30, 100, 10, type = "Var")
#' # Simulates 100 sample of size 30 from a normal distribution with mean 100
#' # and a standard deviation of 10.  From the 100 simulated samples, 95% confidence
#' # intervals for the variance are constructed and depicted in the graph.
#' 
#' cisim(100, 50, 0.5, type = "Pi", conf.level = 0.92)
#' # Simulates 100 samples of size 50 from a binomial distribution where the 
#' # population proportion of successes is 0.5.  From the 100 simulated samples,
#' # 92% confidence intervals for Pi are constructed and depicted in the graph.
#'  
#' @keywords programming 
#####################################################################################
## Fixed back spacing issue and changed matrix approach to a loop approach 6/20/13 ##
#####################################################################################
cisim <- function (samples = 100, n = 30, parameter = 0.5, sigma = 1,
                   conf.level = 0.95, type = c("Mean", "Var", "Pi"))
{
  Adkblue <- "#0080FF"
  Aorange <- "#FF4C0C"
  alpha <- 1 - conf.level
  CL <- conf.level * 100
  N <- samples
  type <- match.arg(type)
  if (length(type) > 1 || is.na(type))
    stop("alternative must be one \"Mean\", \"Var\", \"Pi\"")
  if (type == "Pi" && (parameter <= 0 | parameter >= 1))
    stop("Value for Pi (parameter) must be between 0 and 1.")
  if (N <= 0 || n <= 0)
    stop("Number of random CIs (samples) and sample size (n) must both be at least 1")
  if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) ||
                                 conf.level <= 0 || conf.level >= 1))
    stop("'conf.level' must be a single number between 0 and 1")
  if (sigma <= 0 && (type == "Var" || type == "Mean"))
    stop("Variance must be a positive value")
  ll <- numeric(samples)
  ul <- numeric(samples)
  if (type == "Mean") {
    for(i in 1:samples){
      xbar <- mean(rnorm(n, parameter, sigma))
      ll[i] <- xbar - qnorm(1 - alpha/2)*sigma/sqrt(n)
      ul[i] <- xbar + qnorm(1 - alpha/2)*sigma/sqrt(n)
    }
    notin <- sum((ll > parameter) + (ul < parameter))
    percentage <- round((notin/samples) * 100, 2)
    plot(ll, type = "n", ylim = c(min(ll), max(ul)), xlab = " ",
         ylab = " ")
    title(sub = bquote(paste("Note: ", .(percentage), "% of the random confidence intervals do not contain ", mu, " = ", .(parameter))))
    title(main = bquote(paste(.(samples), " random ", .(CL), "% confidence intervals where ", mu, " = ", .(parameter))))
    for (i in 1:samples) {
      low <- ll[i]
      high <- ul[i]
      if (low < parameter & high > parameter) {
        segments(i, low, i, high)
      }
      else if (low > parameter & high > parameter) {
        segments(i, low, i, high, col = Aorange, lwd = 5)
      }
      else {
        segments(i, low, i, high, col = Adkblue, lwd = 5)
      }
    }
    abline(h = parameter)
    cat(percentage, "\b% of the random confidence intervals do not contain Mu =", parameter, "\b.", "\n")
  }
  else if (type == "Var") {
    for(i in 1:samples){
      s2 <- var(rnorm(n, parameter, sigma))
      ll[i] <- ((n - 1) * s2)/qchisq(1 - alpha/2, (n - 1))
      ul[i] <- ((n - 1) * s2)/qchisq(alpha/2, (n - 1))
    }
    variance <- sigma^2
    notin <- sum((ll > variance) + (ul < variance))
    percentage <- round((notin/samples) * 100, 2)
    plot(ll, type = "n", ylim = c(min(ll), max(ul)), xlab = " ",
         ylab = " ")
    title(sub = bquote(paste("Note: ", .(percentage), "% of the random confidence intervals do not contain ", sigma^2, " = ", .(variance))))
    title(main = bquote(paste(.(samples), " random ", .(CL), "% confidence intervals where ", sigma^2, " = ", .(variance))))
    for (i in 1:samples) {
      low <- ll[i]
      high <- ul[i]
      if (low < variance & high > variance) {
        segments(i, low, i, high)
      }
      else if (low > variance & high > variance) {
        segments(i, low, i, high, col = Aorange, lwd = 5)
      }
      else {
        segments(i, low, i, high, col = Adkblue, lwd = 5)
      }
    }
    abline(h = variance)
    cat(percentage, "\b% of the random confidence intervals do not contain Var =", sigma^2, "\b.", "\n")
  }
  else if (type == "Pi") {
    X <- rbinom(samples, n, parameter)
    p <- X/n
    ll <- p - qnorm(1 - alpha/2) * sqrt((p * (1 - p))/n)
    ul <- p + qnorm(1 - alpha/2) * sqrt((p * (1 - p))/n)
    notin <- sum((ll > parameter) + (ul < parameter))
    percentage <- round((notin/samples) * 100, 2)
    plot(ll, type = "n", ylim = c(min(ll), max(ul)), xlab = " ",
         ylab = " ")
    title(sub = bquote(paste("Note: ", .(percentage), "% of the random confidence intervals do not contain ", pi, " = ", .(parameter))))
    title(main = bquote(paste(.(samples), " random ", .(CL), "% confidence intervals where ", pi, " = ", .(parameter))))
    for (i in 1:samples) {
      low <- ll[i]
      high <- ul[i]
      if (low < parameter & high > parameter) {
        segments(i, low, i, high)
      }
      else if (low > parameter & high > parameter) {
        segments(i, low, i, high, col = Aorange, lwd = 5)
      }
      else {
        segments(i, low, i, high, col = Adkblue, lwd = 5)
      }
    }
    abline(h = parameter)
    cat(percentage, "\b% of the random confidence intervals do not contain Pi =", parameter, "\b.", "\n")
  }
}