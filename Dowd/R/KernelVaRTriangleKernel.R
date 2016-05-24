#' Calculates VaR using triangle kernel approach
#' 
#' The output consists of a scalar VaR for specified confidence level.
#' 
#' @param Ra Profit and Loss data set
#' @param cl VaR confidence level
#' @param plot Bool, plots cdf if true.
#' @return Scalar VaR
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#'
#' @author Dinesh Acharya
#' @examples
#' 
#'    # VaR for specified confidence level using triangle kernel approach
#'    Ra <- rnorm(30)
#'    KernelVaRTriangleKernel(Ra, .95)
#'
#' @export
KernelVaRTriangleKernel <- function(Ra, cl, plot=TRUE) {
  PandL <- as.vector(Ra)
  mu <- mean(PandL)
  sigma <- sd(PandL)
  
  # Obtain pdf values
  kernel.data <- density(PandL, kernel = "triangular", from = mu - 4 * sigma, to = mu + 4 * sigma, n = 1000, bw = "nrd")
  kernel.pdf <- kernel.data$y
  x.values <- kernel.data$x
  delta.x <- x.values[2]-x.values[1]
  n <- 1000 # = length(x.values)
  
  # Obtain cdf values
  cdf <- double(n)
  cdf[1] <- kernel.pdf[1] * delta.x
  for (i in 2:n) {
    cdf[i] <- kernel.pdf[i] * delta.x + cdf[i - 1]
  }
  if(plot == TRUE) {
    plot(x.values, kernel.pdf, type="l", main = "Constructed Pdf")
  }
  
  
  # Derivation of required percentile
  cdf.indices.less.than.prob <- which(cdf<cl)
  # Gives vector of indices for all cdf-values less than probability
  max.cdf.index.less.than.prob <- length(cdf.indices.less.than.prob)
  # Gives index of cdf-value just less than probability
  lower.x.bound <- x.values[max.cdf.index.less.than.prob]
  # Gives x-value just below specified probability
  upper.x.bound <- x.values[max.cdf.index.less.than.prob + 1]
  # Gives x-value just above specified probability
  VaR <- (lower.x.bound + upper.x.bound) / 2 # Desired percentile, ie, answer.
  return(VaR)
  
}