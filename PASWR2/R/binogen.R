#' @title Binomial Distribution Simulation
#' 
#' @description Function that generates and displays \emph{m} repeated samples of \emph{n} Bernoulli trials with a given probability of success
#' 
#' @param samples number of repeated samples to generate
#' @param n number of Bernoulli trials
#' @param pi probability of success for each Bernoulli trial
#' 
#' @return 
#' \item{\code{simulated.distribution}}{Simulated binomial distribution}
#' \item{\code{theoretical.distribution}}{Theoretical binomial distribution}
#' 
#' @author Alan T. Arnholt <arnholtat@@appstate.edu> 
#' 
#' @export
#' 
#' @examples
#' bino.gen(samples=50000, n = 10, pi = 0.80)
#' 
#' @keywords programming
#####################################################################################
bino.gen <- function(samples = 10000, n = 20, pi = 0.5){
  values <- sample(c(0, 1), samples*n, replace=TRUE, prob=c(1 - pi, pi))
  value.mat <- matrix(values, ncol=n)
  Successes <- apply(value.mat, 1, sum)
  a1 <- round((table(Successes)/samples), 3)
  b1 <- round(dbinom(0:n, n, pi), 3)
  names(b1) <- 0:n
  hist(Successes, breaks=c((-.5 + 0):(n + .5)), freq = FALSE, ylab = "", 
       main = " Theoretical Values Superimposed \n Over Histogram of Simulated Values", 
       col = 13, ylim = c(0, max(a1, b1)))
  x <- 0:n
  fx <- dbinom(x, n, pi)
  lines(x, fx, type = "h")
  lines(x, fx, type = "p", pch = 16)
  list(simulated.distribution = a1, theoretical.distribution = b1)
}
