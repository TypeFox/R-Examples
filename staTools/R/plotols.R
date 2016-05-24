#' Plot OLS Fit of Discrete Powerlaw Distributions
#'
#' Plot OLS fit of discrete powerlaw distributions.
#' @param o A discrete powerlaw object.
#' @keywords plot fit powerlaw distributions
#' @import magicaxis
#' @export plotols
#' @examples
#' x = moby
#' o = displo(x)
#' plotols(o)

plotols = function(o)
{
  dev.off()
  d = o$x
  t = as.data.frame(table(d))
  y = log(t$Freq)
  x = log(as.numeric(t$d))
  fit = lm(y ~ x)
  alpha = -unname(fit$coef[2])
  rp = vector('expression',1)
  rp[1] = substitute(expression(hat(alpha) == MYVALUE),
                     list(MYVALUE = format(alpha,dig=3)))[2]
  par(las = 1, mgp = c(1,1,0))
  plot(x,y, col = "gray", pch = 16,
       xlab = "log(x)", ylab = "log(y)", main = "", cex.lab = 1.2,
       xaxt = 'n', yaxt = 'n', xaxs = 'i', yaxs = 'i')
  magaxis(1:2, cex.axis = 1, labels = FALSE)
  abline(fit, col = "blue", lwd = 3)
  legend("bottomleft", rp, lwd = 3, col = "blue", inset = 0.01,
         box.lty = 0, title.adj = 0.1)
}