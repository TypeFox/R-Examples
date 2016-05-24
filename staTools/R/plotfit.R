#' Plot Fit of Discrete Powerlaw Distributions
#'
#' Plot fit of discrete powerlaw distributions.
#' @param o A discrete powerlaw object.
#' @param xmax The maximum value to show.
#' @keywords plot fit powerlaw distributions
#' @import magicaxis
#' @export plotfit
#' @examples
#' x = moby
#' o = displo(x)
#' getXmin(o)
#' plotfit(o)

plotfit = function(o, xmax = 1e5)
{
  dev.off()
  x = o$x
  xmin = o$xmin
  alpha = o$alpha
  obs = o$ux
  ccdf = o$p
  minprob = 1 - pdispl(xmax, xmin = xmin, alpha = alpha)

  rp = vector('expression',2)
  rp[1] = substitute(expression(hat(x)[min] == MYVALUE),
                     list(MYVALUE = format(xmin,dig=0)))[2]
  rp[2] = substitute(expression(hat(alpha) == MYOTHERVALUE),
                     list(MYOTHERVALUE = format(alpha, digits = 3)))[2]

  # plot
  par(las = 1, mgp = c(2.5,1,0))
  plot(obs[-length(obs)],ccdf[-length(ccdf)], log = "xy",
       xlab = "x", ylab = "P(X > x)", main = "",
       cex.lab = 1.2,
       pch = 16, col = "grey",
       xaxt = 'n', yaxt = 'n')
  magaxis(1:2, cex.axis = 1)
  pxmin = 1 - pdispl(xmin, xmin = xmin, alpha = alpha)
  pxmax = 1 - pdispl(xmax, xmin = xmin, alpha = alpha)

  # scale
  d_cdf = o$p #pdis(o$ux, xmin, alpha)
  dif = o$ux - xmin
  upper = which(dif > 0)[1]
  lower = max(upper - 1, 1)
  x_dif = o$ux[lower] - o$ux[upper]
  y_dif = d_cdf[lower] - d_cdf[upper]
  scale = d_cdf[lower] + y_dif*(xmin - o$ux[lower])/x_dif

  lines(c(xmin, xmax),
        c(pxmin*scale, pxmax*scale),
        col = "red", lwd = 3)
  lines(c(obs[(which(obs==xmin)):(length(obs)-1)]),
        c(ccdf[(which(obs==xmin)):(length(obs)-1)]),
        col = "blue", lwd = 3)
  points(xmin,ccdf[which(obs == xmin)], col = "red", pch = 18, cex = 1)
  legend("bottomleft", legend = rp, pch = c(18,NA), col = c("red", "red"),
         lty = c(NA, 1), box.lty = 0, lwd = c(3,3), inset = 0.05)
}
