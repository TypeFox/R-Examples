#' Inspect Discrete Powerlaw Distributions
#'
#' A graphical tool to inspect discrete powerlaw distributions.
#' @param o A discrete powerlaw object.
#' @param plot Logical, whether to show the plot. By default is set to TRUE.
#' @param guess A guess on the true value of the lower bound. By default is set to 1.
#' @param showQ Logical, whether to show the quantiles of the distribution. By default is set to FALSE.
#' @param plothill Logical, whether to show Hill plot. By default is set to TRUE.
#' @param summary Logical, whether to print some information about the powerlaw distribution. By default is set to TRUE.
#' @param xmax The maximum value to consider as candidate for the lower bound.
#' @keywords inspect powerlaw distributions
#' @import magicaxis
#' @export inspect
#' @examples
#' x = moby
#' o = displo(x)
#' inspection = inspect(o, guess = 7)

inspect = function(o, plot = TRUE, guess = 1, showQ = FALSE,
                    plothill = TRUE, summary = TRUE, xmax = 1e5)
{
  inspect = list()
  x = o$x
  nu = o$nux
  obs = o$ux
  ccdf = o$p
  n = o$nx
  quantile = quantile(x, probs = seq(0,1,0.1))
  if (max(x)<xmax) {xmax = max(x)}

  # plot
  if (plot)
  {
    #dev.off()
    nplot = 1
    if (plothill == TRUE){nplot = 2}
    par(las = 1, mgp = c(2.5,1,0), mfrow = c(1,nplot))
    plot(obs[-length(obs)], ccdf[-length(ccdf)], log = "xy",
         xlab = "x", ylab = "P(X > x)", main = "",
         cex.lab = 1.2,
         pch = 16, col = "grey",
         xaxt = 'n', yaxt = 'n')
    if (showQ)
    {
      abline(v = quantile, col = "red", lwd = 2, lty = 3)
    }
    magaxis(1:2, cex.axis = 1)
    if (guess>0) {
      START = obs[which.min(abs(obs - guess))]
      if (START > guess) {
        START = obs[which.min(abs(obs - guess)) - 1]
      }
      guess = obs[which(obs==START)]

      xu = sort(x)
      len_xu = length(xu)
      nq = length(xu[xu>=guess])
      q = xu[(n-nq+1):len_xu]
      q = q[q <= xmax]
      alpha = 1 + length(q) / sum(log(q/(guess-0.5)))

      # scale
      d_cdf = o$p #pdis(o$ux, xmin, alpha)
      dif = o$ux - guess
      upper = which(dif > 0)[1]
      lower = max(upper - 1, 1)
      x_dif = o$ux[lower] - o$ux[upper]
      y_dif = d_cdf[lower] - d_cdf[upper]
      scale = d_cdf[lower] + y_dif*(guess - o$ux[lower])/x_dif

      pxmin = 1 - pdispl(guess, xmin = guess, alpha = alpha)
      pxmax = 1 - pdispl(xmax, xmin = guess, alpha = alpha)
      lines(c(guess, xmax),
            c(pxmin*scale, pxmax*scale),
            col = "red", lwd = 3)
      #       lines(c(guess, tail(obs,2)[1]),
      #             c(ccdf[which(obs == guess)], ccdf[which(obs == tail(obs,2)[1])]),
      #             col = "blue", lwd = 3)
      lines(c(obs[(which(obs==guess)):(length(obs)-1)]),
            c(ccdf[(which(obs==guess)):(length(obs)-1)]),
            col = "blue", lwd = 3)
      points(guess,ccdf[which(obs == guess)], col = "red", pch = 18, cex = 1)

      if (plothill == TRUE)
        plothill(o, gxmin = guess, xmax = xmax)
    }
  }

  if (summary)
  {
    cat("\nSUMMARY",
        "\n*******",
        "\nN = ", n,
        "\n(unique) = ", nu,
        "\nmax value = ", max(x),
        "\nxmin (guess) = ", guess,
        "\nalpha (guess) = ", alpha)
  }

  # to return
  inspect$n = n
  inspect$nu = nu
  inspect$obs = obs
  inspect$prob = ccdf
  inspect$quantile = quantile
  inspect$xmax = xmax
  inspect$gxmin = guess
  inspect$galpha = alpha
  return(inspect)
}
