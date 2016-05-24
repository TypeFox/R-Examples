#' @title Simulated Distribution of \eqn{D_n}{D[n]} (Kolmogorov-Smirnov)
#' 
#' @description Function to visualize the sampling distribution of \eqn{D_n}{D[n]} (the Kolmogorov-Smirnov one sample statistic) and to find simulated critical values.
#' 
#' @param n sample size
#' @param sims number of simulations to perform
#' @param alpha desired \eqn{\alpha} level
#' 
#' @author Alan T. Arnholt <arnholtat@@appstate.edu> 
#'  
#' @seealso \code{\link{ksldist}}
#' @export
#' 
#' @examples
#' ksdist(n = 10, sims = 15000, alpha =0.05)
#' @keywords hplot
###########################################################################
ksdist <- function (n = 10, sims = 10000, alpha = 0.05)
{
  Dn <- replicate(sims, ks.test(rnorm(n), pnorm)$statistic)
  cv <- quantile(Dn, 1 - alpha)
  plot(density(Dn), col = "blue", lwd = 2, main = "",
       xlab = paste("Simulated critical value =", round(cv,3),
                    "for n =", n, "when the alpha value =", alpha))
  title(main = list(expression(paste("Simulated Sampling Distribution of " ,
                                     D[n]))))
}