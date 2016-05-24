#' Histogram and density plot
#' 
#' This function plots a histogram and a density plot of a single variable 
#' using the R-package \pkg{lattice}.
#' 
#' @param x Vector of numeric values of the variable for which the histogram 
#' nd density plot should be created.
#' 
#' @param HD Character value indicating the type of plot to be created. 
#' Available options are \code{"over"}, to create a histogram superimposed by 
#' the theoretical density plot of a normally distributed variable, and 
#' \code{"stack"}, to create a histogram and an empirical density plot in
#' separated panels. Defaults to \code{HD = "over"}.
#' 
#' @param nint Integer specifying the number of histogram bins. Defaults to
#' \code{nint = 20}.
#' 
#' @param stats Logical to indicate if descriptive statistics of the variable
#' \code{x} should be added to the plot. Available only when \code{HD = "over"}.
#' The function tries to automatically find the best location to put the 
#' descriptive statistics given the shape of the histogram. Defaults to 
#' \code{stats = TRUE}.
#' 
#' @param digits Integer indicating the number of decimal places to be used 
#' when printing the statistics of the variable \code{x}. Defaults to 
#' \code{digits = 2}.
#' 
#' @param BoxCox Logical to indicate if the variable \code{x} should be 
#' transformed using the Box-Cox family of power transformations. The estimated
#' lambda value of the Box-Cox transform is printed in the console. It is set to
#' zero when negative. Defaults to \code{BoxCox = FALSE}.
#' 
#' @param col Vector of two elements, the first indicating the colour of the
#' histogram, the second indicating the colour of the density plot. Defaults to 
#' \code{col = c("lightgray", "black")}.
#' 
#' @param lwd Vector of two elements, the first indicating the line width of the
#' histogram, the second indicating the line width of the density plot. Defaults 
#' to \code{lwd = c(1, 1)}.
#' 
#' @param lty Character value indicating the line type for the density plot.
#' Defaults to \code{lty = "dashed"}.
#' 
#' @param xlim Vector of two elements defining the limits of the x axis. The
#' function automatically optimizes \code{xlim} based on the density plot.
#' 
#' @param ylim Vector of two elements defining the limits of the y axis. The
#' function automatically optimizes \code{ylim} based both histogram and density 
#' plot.
#' 
#' @param ... Other arguments that can be passed to \pkg{lattice} functions. 
#' There is no guarantee that they will work.
#' 
#' @details
#' The user should visit the help pages of \code{\link[lattice]{histogram}},
#' \code{\link[lattice]{densityplot}}, \code{\link[lattice]{panel.mathdensity}}, 
#' \code{\link[car]{powerTransform}} and \code{\link[car]{bcPower}} to obtain 
#' more details about the main functions used to built \code{plotHD}.
#' 
#' @return
#' An object of class \code{"trellis"}. The 
#' \code{\link[lattice]{update.trellis}} method can be used to update components
#' of the object and the \code{\link[lattice]{print.trellis}} print method 
#' (usually called by default) will plot it on an appropriate plotting device.
#' 
#' @references
#' Sarkar, Deepayan (2008) \emph{Lattice: Multivariate Data Visualization with 
#' R}, Springer. \url{http://lmdvr.r-forge.r-project.org/}
#' 
#' @author Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' 
#' @seealso \code{\link[lattice]{histogram}}, 
#' \code{\link[lattice]{densityplot}}, \code{\link[lattice]{panel.mathdensity}},
#' \code{\link[car]{powerTransform}}, \code{\link[car]{bcPower}}.
#' @import lattice latticeExtra
#' @importFrom stats update
#' @export
#' @examples
#' x <- rnorm(100, 10, 2)
#' plotHD(x, HD = "stack")
#' plotHD(x, HD = "over")
#' 
#' @keywords dplot
#' 
# FUNCTION #####################################################################
plotHD <- 
  function (x, HD = "over", nint = 20, digits = 2, stats = TRUE, 
            BoxCox = FALSE, col = c("lightgray", "black"), lwd = c(1, 1),
            lty = "dashed", xlim, ylim, ...) {
    
    # Check if suggested packages are installed
    if (!requireNamespace("moments", quietly = TRUE)) {
      stop(paste("Package 'moments' needed for this function to work. ",
                 "Please install it.", sep = ""), call. = FALSE)
    }
    if (!requireNamespace("car", quietly = TRUE)) {
      stop(paste("Package 'car' needed for this function to work. ",
                 "Please install it.", sep = ""), call. = FALSE)
    }
    
    if (BoxCox) {
      
      # Check if the variable has negative values
      check <- any(x <= 0)
      if (check) {
        message("data has negative values...")
        x <- x + abs(min(x)) + 1
      }
      
      lambda <- car::powerTransform(x)
      print(summary(lambda))
      lambda <- as.numeric(lambda$lambda)
      if (lambda < 0) {
        message("estimated lambda value is negative... setting to zero")
        lambda <- 0
      }
      x <- car::bcPower(x, lambda)
    } else {
      lambda <- 1
    }
    if (HD == "over") {
      #if (missing(xlim)) {
      #  xlim <- densityplot(x)$x.limits
      #}
      p <- lattice::histogram(
        x, type = "density", col = col[1], xlim = xlim, nint = nint, 
        lwd = lwd[1], ..., 
        panel = function(x, ...) {
          lattice::panel.histogram(x, ...)
          lattice::panel.rug(x, col = col[2], lwd = lwd[1])
          lattice::panel.mathdensity(dmath = stats::dnorm, col = col[2],
                                     lwd = lwd[2], lty = lty, 
                                     args = list(mean = mean(x), 
                                                 sd = stats::sd(x)),
                                     n = length(x))
          })
      if (missing(ylim)) {
        y1 <- p$y.limits
        y2 <- lattice::densityplot(x)$y.limits
        ylim <- c(min(c(y1[1], y2[1])), max(c(y1[2], y2[2])))
      }
      p$y.limits <- ylim
      if (stats) {
        skw <- round(c(moments::skewness(x)), 2)
        leg <- c(paste("Lambda = ", round(lambda, 4), "\n", 
                       "Mean = ", round(mean(x), digits), " (", 
                                  round(stats::sd(x), digits), ")\n",
                       "Median = ", round(stats::median(x), digits), "\n",
                       "Range = ", round(min(x), digits), "-", 
                                   round(max(x), digits), "\n",
                       "Skew = ", skw, sep = ""))
        if (skw >= 1) {
          y <- NA
          pos <- NA
          p <- p + 
            latticeExtra::layer(lattice::panel.text(x = x, y = y, labels = leg,
                                                    pos = pos),
                                data = list(x = c(max(p$x.limits) * 0.99), 
                                            y = c(max(p$y.limits) * 0.9), 
                                            leg = leg, pos = 2))
        }
        if (skw < 1) {
          y <- NA
          pos <- NA
          p <- p + 
            latticeExtra::layer(lattice::panel.text(x = x, y = y, labels = leg, 
                                                    pos = pos),
                                data = list(x = c(min(p$x.limits)),
                                            y = c(max(p$y.limits) * 0.9),
                                            leg = leg, pos = 4))
        }
      }
    }
    if (HD == "stack") {
      p2 <- lattice::densityplot(x, col = col[2], pch = 20, cex = 0.5, 
                                 n = length(x), lwd = lwd[2], lty = lty, ...)
      if (missing(xlim)) {
        xlim <- p2$x.limits
      }
      p1 <- lattice::histogram(x, col = col[1], xlim = xlim, nint = nint, 
                               lwd = lwd[1], ...)
      p <- update(c(p1, p2), layout = c(1, 2), 
                  ylab = list(c("Percent of Total", "Density")), ...)
    }
    return(p)
  }
# End!
