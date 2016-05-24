#' @title Normal Area
#' 
#' @description Function that computes and draws the area between two user specified values in a user specified normal distribution with a given mean and standard deviation
#' 
#' @param lower the desired lower value 
#' @param upper the desired upper value
#' @param m the mean for the population (default is the standard normal with \code{m = 0})
#' @param sig the standard deviation of the population (default is the standard normal with \code{sig = 1})
#' 
#' @return Draws the specified area in a graphics device
#' 
#' @author Alan T. Arnholt <arnholtat@@appstate.edu> 
#' 
#' @export
#' 
#' @examples
#' # Finds and graphically illustrates P(70 < X < 130) given X is N(100, 15)
#' normarea(lower = 70, upper = 130, m = 100, sig = 15) 
#' 
#' @keywords hplot
########################################################################
normarea <- function (lower = -Inf, upper = Inf, m = 0, sig = 1)
{
  Altblue <- "#CDCDED"
  Fontcol <- "#3333B3"
  opar <- par(no.readonly = TRUE)
  par(mar=c(4, 1, 4, 1))
  area <- pnorm(upper, m, sig) - pnorm(lower, m, sig)
  ra <- round(area, 4)
  x <- seq(m - 4 * sig, m + 4 * sig, length = 1000)
  y <- dnorm(x, m, sig)
  par(pty = "m")
  plot(x, y, type = "n", xaxt = "n", yaxt = "n", xlab = "",
       ylab = "", main = "")
  
  mtext(substitute(paste("The area between ",lower," and ",upper," is ",ra)),
        side = 3, line = 1, font = 2, cex = 1.15)
  mtext(substitute(paste("RV ~ N (" ,mu == m,", ",sigma == sig,")" ),
                   list(m = m, sig = sig)), side = 1, line = 3, col = Fontcol)
  
  if (lower == -Inf || lower < m - 4 * sig) {
    lower <- m - 4 * sig
  }
  if (upper == Inf || upper > m + 4 * sig) {
    upper <- m + 4 * sig
  }
  axis(1, at = c(m, lower, upper), labels = c(m, lower, upper))
  xaxis1 <- seq(lower, upper, length = 200)
  yaxis1 <- dnorm(xaxis1, m, sig)
  xaxis1 <- c(lower, xaxis1, upper)
  yaxis1 <- c(0, yaxis1, 0)
  polygon(xaxis1, yaxis1, density = -1, col = Altblue)
  lines(x, y, lwd = 2)
  lines(c(m - 4 * sig, m + 4 * sig), c(0, 0), lwd = 2)
  lines(c(lower, lower), c(0, dnorm(lower, m, sig)), lwd = 2)
  lines(c(upper, upper), c(0, dnorm(upper, m, sig)), lwd = 2)
  on.exit(par(opar))
}