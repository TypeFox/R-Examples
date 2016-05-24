#' @title Interval Plot
#' 
#' @description Function to graph intervals
#' 
#' @param ll vector of lower values
#' @param ul vector of upper values
#' @param parameter value of the desired parameter (used when graphing confidence intervals)
#' 
#' @return Draws user-given intervals on a graphical device.
#' 
#' @author Alan T. Arnholt <arnholtat@@appstate.edu> 
#' 
#' @export
#' 
#' @examples
#' set.seed(385)
#' samples <- 100
#' n <- 625
#' ll <- numeric(samples)
#' ul <- numeric(samples)
#' xbar <- numeric(samples)
#' for (i in 1:samples){
#'   xbar[i] <- mean(rnorm(n, 80, 25))
#'   ll[i] <- xbar[i] - qnorm(.975)*25/sqrt(n)
#'   ul[i] <- xbar[i] + qnorm(.975)*25/sqrt(n)
#'   }
#' interval.plot(ll, ul, parameter = 80)
#' 
#' @keywords programming 
#' 
#############################################################################
# interval.plot() recoded to be more general...6/27/13
#############################################################################
interval.plot <- function(ll, ul, parameter = 0){
  Adkblue <- "#0080FF"
  Aorange <- "#FF4C0C"
  N <- length(ll)
  notin <- sum((ll > parameter) + (ul < parameter))
  percentage <- round((notin/N) * 100, 2)
  plot(ll, type = "n", ylim = c(min(ll), max(ul)), xlab = " ", ylab = " ")
  title(sub = bquote(paste("Note: ", .(percentage), 
                           "% of the intervals do not contain the parameter = ", 
                           .(parameter))))
  for (i in 1:N) {
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
  cat(percentage, "\b% of the intervals do not contain the parameter =", 
      parameter, "\b.", "\n")
}

