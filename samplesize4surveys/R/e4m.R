#' @import TeachingSampling
#' @export
#' 
#' @title
#' Statistical errors for the estimation of a single mean 
#' @description 
#' This function computes the cofficient of variation and the standard error when estimating a single mean under a complex sample design.
#' @return 
#' The coefficient of variation and the margin of error for a predefined sample size.
#' @details
#' We note that the coefficent of variation is defined as: \deqn{cve = \frac{\sqrt{Var(\bar{y}_S)}}{\bar{y}_S}} 
#' Also, note that the magin of error is defined as: \deqn{\varepsilon = z_{1-\frac{\alpha}{2}}\sqrt{Var(\bar{y}_S)}}
#' 
#' @author Hugo Andres Gutierrez Rojas <hugogutierrez at usantotomas.edu.co>
#' @param N The population size.
#' @param n The sample size.
#' @param mu The value of the estimated mean of the variable of interest.
#' @param sigma The value of the standard deviation of the variable of interest.
#' @param DEFF The design effect of the sample design. By default \code{DEFF = 1}, which corresponds to a simple random sampling design.
#' @param conf The statistical confidence. By default \code{conf = 0.95}.
#' @param plot Optionally plot the errors (cve and margin of error) against the sample size.
#' 
#' @references 
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas y estimacion de parametros}. Editorial Universidad Santo Tomas
#' @seealso \code{\link{ss4p}}
#' @examples 
#' e4m(N=10000, n=400, mu = 10, sigma = 10)
#' e4m(N=10000, n=400, mu = 10, sigma = 10, plot=TRUE)
#' e4m(N=10000, n=400, mu = 10, sigma = 10, DEFF=3.45, conf=0.99, plot=TRUE)


e4m <- function(N, n, mu, sigma, DEFF = 1, conf = 0.95, plot = FALSE) {
  
  S2 <- sigma^2 * DEFF
  Z <- 1 - ((1 - conf)/2)
  f <- n/N
  VAR <- (1/n) * (1 - f) * sigma^2
  CVE <- 100 * sqrt(VAR)/mu
  ME <- qnorm(Z) * sqrt(VAR)
  
  if (plot == TRUE) {
    
    nseq <- seq(1, N, 10)
    cveseq <- rep(NA, length(nseq))
    meseq <- rep(NA, length(nseq))
    
    for (k in 1:length(nseq)) {
      fseq <- nseq[k]/N
      varseq <- (1/nseq[k]) * (1 - fseq) * S2
      cveseq[k] <- 100 * sqrt(varseq)/mu
      meseq[k] <- qnorm(Z) * sqrt(varseq)
    }
    
    par(mfrow = c(1, 2))
    plot(nseq, cveseq, type = "l", lty = 1, pch = 1, col = 3, ylab = "Coefficient of variation (%)", xlab = "Sample Size")
    points(n, CVE, pch = 8, bg = "blue")
    abline(h = CVE, lty = 3)
    abline(v = n, lty = 3)
    
    plot(nseq, meseq, type = "l", lty = 1, pch = 1, col = 3, ylab = "Margin of error", xlab = "Sample Size")
    points(n, ME, pch = 8, bg = "blue")
    abline(h = ME, lty = 3)
    abline(v = n, lty = 3)
  }
  
  msg <- cat("With the parameters of this function: N =", N, "n = ", n, "mu =", mu,  "sigma =", sigma,
             "DEFF = ", DEFF, "conf =", conf, ". \nThe estimated coefficient of variation is ", 
             CVE, ". \nThe margin of error is", ME, ". \n \n")
  
  result <- list(cve = CVE, Margin_of_error = ME)
  result
} 