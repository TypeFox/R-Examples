###########################################################################
##                                                                       ##
## plot.dissimilarities - 'plot' method for class 'dissimilarities'      ##
##                                                                       ##
## Created       : 27-May-2006                                           ##
## Author        : Gavin Simpson                                         ##
## Version       : 0.1                                                   ##
## Last modified : 27-May-2006                                           ##
##                                                                       ##
## ARGUMENTS:                                                            ##
## x             - object on which method dispatch applied               ##
## prob          - density probability for close modern analogue         ##
## legend        - logical. Draw a legend                                ##
## n.rnorm       - number of random normal deviates for reference line   ##
## col, col.ref  - colours for the dissimilarity and reference density   ##
##                 functions drawn                                       ##
## lty, lty.ref  - line types for the dissimilarity and reference        ##
##                 density functions drawn                               ##
## xlab, ylab    - x- and y-axis labels                                  ##
## main, sub     - main and subtitle for the plot                        ##
## ...           - arguments passed to other graphics functions          ##
##                                                                       ##
###########################################################################
plot.dissimilarities <- function(x,
                                 prob = 0.05,
                                 legend = TRUE,
                                 n.rnorm = 100000,
                                 col = "black",
                                 col.ref = "red",
                                 lty = "solid",
                                 lty.quant = "dotted",
                                 xlab = NULL,
                                 ylab = NULL,
                                 main = NULL,
                                 sub = NULL, ...)
  {
    xlab <- if(is.null(xlab))
      "Dissimilarity"
    ylab <- if(is.null(ylab))
      "Density"
    main <- if(is.null(main))
      "Distribution of training set dissimilarities"
    distr.mean <- mean(x)
    distr.sd <- sd(x)
    len <- length(x)
    distr.dens <- density(x)
    ref <- rnorm(n.rnorm, mean = distr.mean, sd = distr.sd)
    ref.dens <- density(ref)
    sub <- if(is.null(sub))
      paste("Observed N =", distr.dens$n, " Bandwidth =", round(distr.dens$bw, 3))
    ylims <- range(distr.dens$y, ref.dens$y)
    xlims <- range(distr.dens$x, ref.dens$x)
    plot(distr.dens, xlim = xlims, ylim = ylims, xlab = xlab, col = col,
         sub = sub, main = main)
    lines(ref.dens, col = col.ref)
    distr.quan <- quantile(x, prob = prob)
    ref.quan <- quantile(ref, prob = prob)
    abline(v = distr.quan, lty = lty.quant[1], col = col)
    abline(v = ref.quan, lty = lty.quant[1], col = col.ref)
    if(legend)
      {
        legend("topright", inset = 0.02,
               legend = c("Observed", "Reference",
                 paste(prob*100, "th % (Obs. == ",
                       round(distr.quan, 3),  ")", sep = ""),
                 paste(prob*100, "th % (Ref. == ",
                       round(ref.quan,3), ")", sep = "")),
               col = rep(c(col, col.ref), 2),
               lty = rep(c(lty, lty.quant), each = 2),
               cex = 0.6
               )
      }
    invisible()
  }
