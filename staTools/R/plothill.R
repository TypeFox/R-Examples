#' Hill plot
#'
#' Hill plot for discrete power law distributions.
#' @param o A discrete powerlaw object.
#' @param gxmin Guess on the true value of the lower bound.
#' @param xmax Maximum value considered as candidate for the lower bound. Default is set to 1e5.
#' @keywords plot hill
#' @import magicaxis
#' @export plothill
#' @examples
#' x = moby
#' o = displo(x)
#' plothill(o)

plothill = function(o, gxmin = 0, xmax = 1e5)
{
  if (length(o$xmins) == 0 || length(o$alphas) == 0) {
  x = o$x
  xu = x
  len_xu = length(xu)
  xmin = o$ux
  xmin = xmin[xmin <= xmax]
  alpha = rep(0,(length(xmin)-1))
  n = length(x)
  for (i in 1:(length(xmin)-1))
  {
    nq = length(xu[xu>=xmin[i]])
    q = xu[(n-nq+1):len_xu]
    q = q[q <= xmax]
    alpha[i] = 1 + length(q) / sum(log(q/(xmin[i]-0.5)))
  }

  o$xmins = xmin
  o$alphas = alpha
} else {
  xmin = o$xmins
  alpha = o$alphas
}
  galpha = alpha[which(xmin == gxmin)]
  par(las = 1, mgp = c(2.25,1,0))
  plot(xmin[-length(xmin)],alpha, log = "x", ylim = c(1,3),
       cex.lab = 1.2, xaxt = 'n', yaxt = 'n', xaxs = 'i', yaxs = 'i',
       xlab = expression(x[min]), ylab = expression(alpha),
       pch = 16, col = "grey", cex = 0.75)
  if (gxmin > 0 && galpha > 0) {
    abline(v = gxmin, lwd = 1, col = "red", lty = 2)
    abline(h = galpha, lwd = 1, col = "red", lty = 2)
    points(gxmin, galpha, pch = 16, col = "red", cex = 0.75)
    points(gxmin, galpha, pch = 1, col = "red", cex = 2)
  }
  magaxis(1:2, cex.axis = 1)
}
