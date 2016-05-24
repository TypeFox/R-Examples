#' Plot Multiple Fit of Discrete Powerlaw Distributions
#'
#' Plot multiple fit of discrete powerlaw distributions.
#' @param o A discrete powerlaw object.
#' @keywords plot fit powerlaw distributions
#' @import magicaxis
#' @export plotmultifit
#' @examples
#' x = moby
#' o = displo(x)
#' getXmin(o)
#' plotmultifit(o)

plotmultifit = function(o)
{
  x = o$x
  xmin = o$xmin
  alpha = o$alpha
  obs = o$ux
  ccdf = o$p
  xmax = o$xmax

  x = x[x>xmin]
  # fit exponential
  o$fitexp = list()
  o$fitexp[["rate"]] = 1/mean(x)

  # fit poisson
  o$fitpois = list()
  o$fitpois[["lambda"]] = mean(x)

  # fit lognormal
  o$fitlnorm = list()
  o$fitlnorm[["meanlog"]] = sum(log(x))/ length(x)
  o$fitlnorm[["sdlog"]] = sqrt(sum((log(x) - o$fitlnorm[["meanlog"]])^2)/length(x))

  rp = vector('expression', 5)
  rp[1] = "powerlaw fit"
  rp[2] = "exponential fit"
  rp[3] = "lognormal fit"
  rp[4] = "poisson fit"
  rp[5] = "empirical fit"

  # plot
  par(las = 1, mgp = c(2.5,1,0))
  plot(obs[-length(obs)],ccdf[-length(ccdf)], log = "xy",
       xlab = "x", ylab = "P(X > x)", main = "",
       cex.lab = 1.2,
       pch = 16, col = "grey",
       xaxt = 'n', yaxt = 'n')
  magaxis(1:2, cex.axis = 1)

  # powerlaw fit
  d_cdf = o$p
  dif = o$ux - xmin
  upper = which(dif > 0)[1]
  lower = max(upper - 1, 1)
  x_dif = o$ux[lower] - o$ux[upper]
  y_dif = d_cdf[lower] - d_cdf[upper]
  scale = d_cdf[lower] + y_dif*(xmin - o$ux[lower])/x_dif

  pxmin = 1 - pdispl(xmin, xmin = xmin, alpha = alpha)
  pxmax = 1 - pdispl(xmax, xmin = xmin, alpha = alpha)
  lines(c(xmin, xmax),
        c(pxmin*scale, pxmax*scale),
        col = "red", lwd = 3)
  lines(c(obs[(which(obs==xmin)):(length(obs)-1)]),
        c(ccdf[(which(obs==xmin)):(length(obs)-1)]),
        col = "blue", lwd = 3)
  points(xmin,ccdf[which(obs == xmin)], col = "red", pch = 18, cex = 1)

  # exponential fit
  pxmin = 1 - pexp(min(x), rate = o$fitexp$rate)
  xmax = x[which((1 - pexp(x, o$fitexp$rate)) > 0)]
  xmax = xmax[length(xmax)]
  pxmax = 1 - pexp(xmax, o$fitexp$rate)
  lines(c(xmin, xmax),
        c(pxmin*scale, pxmax*scale),
        col = "purple", lwd = 3)

  # lognormal fit
  pxmin = 1 - plnorm(min(x), meanlog = o$fitlnorm$meanlog, sdlog = o$fitlnorm$sdlog)
  xmax = x[which((1 - plnorm(x, meanlog = o$fitlnorm$meanlog, sdlog = o$fitlnorm$sdlog)) > 0)]
  xmax = xmax[length(xmax)]
  pxmax = 1 - plnorm(xmax, meanlog = o$fitlnorm$meanlog, sdlog = o$fitlnorm$sdlog)
  lines(c(xmin, xmax),
        c(pxmin*scale, pxmax*scale),
        col = "black", lwd = 3)

  # poisson fit
  pxmin = 1 - ppois(min(x), lambda = o$fitpois$lambda)
  xmax = x[which((1 - ppois(x, lambda = o$fitpois$lambda)) > 0)]
  xmax = xmax[length(xmax)]
  pxmax = 1 - ppois(xmax, lambda = o$fitpois$lambda)
  lines(c(xmin, xmax),
        c(pxmin*scale, pxmax*scale),
        col = "forestgreen", lwd = 3)

  legend("bottomleft", legend = rp,
         col = c("red", "purple","black","forestgreen","blue"),
         pch = c(NA,NA,NA,NA,NA),
         lty = c(1,1,1,1,1), box.lty = 0,
         lwd = 3, inset = 0.01, cex = 0.75)
}