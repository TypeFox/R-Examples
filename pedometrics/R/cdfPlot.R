#' Plot estimated cumulative distribution function with confidence limits
#' 
#' This function is a modified version of \code{cdf.plot()} of
#' \pkg{spsurvey}-package including new argument options.
#' 
#' Parameter \code{type.plot} is used only when \code{type.cdf = "Continuous"}.
#' 
#' Care should be taken with possible conflicts between the arguments of the
#' original function \code{\link[spsurvey]{cdf.plot}} and those passed to
#' \code{plot()} using \code{...}. The existence of conflicts between these two
#' functions was one of the reasons for creating this new implementation.
#' 
#' @param units.cdf Indicator for the type of units in which the CDF is
#' plotted, where \dQuote{percent} means the plot is in terms of percent of the
#' population, and \dQuote{units} means the plot is in terms of units of the
#' population.  Defaults to \code{units.cdf = "percent"}.
#' @param type.cdf Character string consisting of the value \dQuote{continuous}
#' or \dQuote{ordinal} that controls the type of CDF plot for each indicator.
#' Defaults to \code{type.cdf = "continuous"}.
#' @param logx Character string consisting of the value \code{""} or \code{"x"}
#' that controls whether the x axis uses the original scale (\code{""}) or the
#' base 10 logarithmic scale (\code{"x"}). Defaults to \code{logx = ""}.
#' @param xlbl Character string providing the x-axis label. If this argument
#' equals \code{NULL}, then the indicator name is used as the label. Defaults
#' to \code{xlbl = NULL}.
#' @param ylbl Character string providing the the y-axis label. Defaults to
#' \code{ylbl = "Percent"}.
#' @param ylbl.r Character string providing the label for the right side
#' y-axis, where \code{ylbl.r = NULL} means a label is not created, and
#' \code{ylbl.r = "Same"} means the label is the same as the left side label
#' (i.e., argument \code{ylbl}). Defaults to \code{ylbl.r = NULL}.
#' @param figlab Character string providing the plot title. Defaults to
#' \code{figlab = NULL}.
#' @param legloc Indicator for location of the plot legend, where \code{legloc
#' = "BR"} means bottom right, \code{legloc = "BL"} means bottom left,
#' \code{legloc = "TR"} means top right, and \code{legloc = "TL"} means top
#' left. Defaults to \code{legloc = "BR"}.
#' @param confcut Numeric value that controls plotting confidence limits at the
#' CDF extremes. Confidence limits for CDF values (percent scale) less than
#' \code{confcut} or greater than 100 minus \code{confcut} are not plotted.  A
#' value of zero means confidence limits are plotted for the complete range of
#' the CDF. Defaults to \code{confcut = 5}.
#' @param conflev Numeric value of the confidence level used for confidence
#' limits. Defaults to \code{conflev = 95}.
#' @param ...  Additional arguments passed to \code{plot()}. See
#' \sQuote{Details}.
#' @param obj Object with the estimated CDF. The resulting object of
#' \code{cont.analysis()} of \pkg{spsurvey}-package.
#' @param ind Indicator variable. The name of the variable as displayed in the
#' resulting object of \code{cont.analysis()}.
#' @param type.plot Type of plot. Desired type of plot to be produced, with
#' options \code{type.plot = "l"}, for \sQuote{line}, and \code{type.plot =
#' "s"} for \sQuote{stair}. See \sQuote{Details}. Defaults to \code{type.plot =
#' "s"}.
#' @param show.param Logical for showing the parameters of the CDF. Available
#' parameters are the mean, the median, and a percentile defined by the
#' argument \code{conflev}.  The legend displays de actual values of all three
#' parameters, including the standard deviation of the mean. The percentile
#' value is calculated using \code{spsurvey::interp.cdf()}.
#' @param round Numeric to set the rounding level of the parameters of the CDF.
#' @param col.param Color of the lines showing the parameters of the CDF.
#' Defaults to \code{col.param = "black"}.
#' @param show.conflev Logical for showing the confidence limits of the CDF.
#' Defaults to \code{show.conflev = TRUE}.
#' @return A plot of the estimated cumulative distribution function with
#' confidence limits.
#' @note Most of the source code that constitutes this function was originaly
#' published in the \pkg{spsurvey}-package, version 2.6 (2013-09-20). The
#' authors were asked to include a few new functionalities, but did not seem to
#' be interested in doing so, since no reply was obtained. This implementation
#' is a way of including such functionalities. When using this function, credit
#' should be given to the authors of the original implementation in the
#' \pkg{spsurvey}-package.
#' @author Tony Olsen \email{Olsen.Tony@@epa.gov}\cr Tom Kincaid
#' \email{Kincaid.Tom@@epa.gov}\cr Alessandro Samuel-Rosa
#' \email{alessandrosamuelrosa@@gmail.com}
#' @seealso \code{\link[spsurvey]{cdf.plot}}.
#' @references Brus, D. J., Kempen, B. and Heuvelink, G. B. M. (2011).
#' Sampling for validation of digital soil maps.  \emph{European Journal of
#' Soil Science}, v. 62, p. 394-407.
#' 
#' Diaz-Ramos, S., D.L. Stevens, Jr., and A.R. Olsen. (1996).  \emph{EMAP
#' Statistical Methods Manual}. EPA/620/R-96/XXX. Corvallis, OR: U.S.
#' Environmental Protection Agency, Office of Research and Development,
#' National Health Effects and Environmental Research Laboratory, Western
#' Ecology Division.
#' 
#' Kincaid, T. M. and Olsen, A. R. (2013) \emph{spsurvey: Spatial Survey Design
#' and Analysis}.  R package version 2.6. URL:
#' \url{http://www.epa.gov/nheerl/arm/}.
#' @keywords dplot hplot
#' @export
#' @examples
#' 
#' \dontrun{
#' ## Estimate the CDF
#' my.cdf <- spsurvey::cont.analysis(spsurvey.obj = my.spsurvey)
#' 
#' ## See indicator levels in the resulting object
#' levels(my.cdf$Pct$Indicator)
#' 
#' ## Plot CDF
#' cdfPlot(obj = my.cdf, ind = "dz", figlab = "",
#'    xlbl = "Difference (m)", xlim = c(-30, 10), type.plot = "s")
#' }
#' 
# FUNCTION #####################################################################
cdfPlot <- 
  function (obj, ind, units.cdf = "percent", type.plot = "s", 
            type.cdf = "continuous", logx = "", xlbl = NULL, ylbl = "Percent",
            ylbl.r = NULL, figlab = NULL, legloc = "BR", confcut = 5,
            show.conflev = TRUE,
            conflev = 95, show.param = TRUE, round = 0, 
            col.param = "black", ...) {
    
    # Check if suggested packages are installed
    if (!requireNamespace("spsurvey", quietly = TRUE)) {
      stop(paste("Package 'spsurvey' needed for this function to work. ",
                 "Please install it.", sep = ""), call. = FALSE)
    }
    
    op <- graphics::par(mgp = c(1.7, 0.6, 0), mar = c(3, 3, 2, 4) + 0.1)
    obj <- obj
    ind <- ind
    cdfest <- obj$CDF[obj$CDF$Indicator == paste(ind), ]
    if (units.cdf == "percent") {
      cdfdata <- cdfest[, c(4, 6, 8, 9, 10)]
    } else if (units.cdf == "units") {
      cdfdata <- cdfest[, c(4, 10, 12, 13, 6)]
    } else {
      stop(paste("the choice of units must be either 'percent' or 'units'"))
    }
    pctval <- c(confcut, 100 - confcut)
    tvalue <- cdfest[, 6] >= pctval[1] & cdfest[, 6] <= pctval[2]
    x <- spsurvey::interp.cdf(pctval, cdfest[, 6], cdfdata[, 1])
    ylow <- spsurvey::interp.cdf(pctval, cdfest[, 6], cdfdata[, 3])
    yhi <- spsurvey::interp.cdf(pctval, cdfest[, 6], cdfdata[, 4])
    if (units.cdf == "percent") {
      ylimit <- c(0, 100)
    } else if (units.cdf == "units") {
      ylimit <- pretty(c(min(c(cdfdata[, 2], ylow)), max(c(cdfdata[, 2], yhi))))
      ylimit <- ylimit[c(1, length(ylimit))]
    }
    if (type.cdf == "continuous") {
      ty <- c("l", "s")
      if(any(ty == type.plot)){
        graphics::plot(cdfdata[, 1], cdfdata[, 2], type = type.plot, 
                       ylim = ylimit, xlab = xlbl, ylab = ylbl, log = logx, ...)
        value <- c(x[1], cdfdata[, 1][tvalue], x[2])
        lower <- c(ylow[1], cdfdata[, 3][tvalue], ylow[2])
        upper <- c(yhi[1], cdfdata[, 4][tvalue], yhi[2])
        if (show.conflev) {  # Logical for showing the CDF's confidence limits
          graphics::lines(value, lower, type = type.plot, lty = 3, lwd = 1.5)
          graphics::lines(value, upper, type = type.plot, lty = 3, lwd = 1.5)
        }
        if (show.param) {
          mea <- round(mean(cdfdata[, 1]), round)
          med <- round(stats::median(cdfdata[, 1]), round)
          per <- round(spsurvey::interp.cdf(conflev, cdfest[, 6], cdfdata[, 1]), round)
          graphics::lines(rep(mea, 2), ylimit, lty = "dashed", col = col.param)
          graphics::lines(rep(med, 2), ylimit, lty = "dotdash", col = col.param)
          graphics::lines(rep(per, 2), ylimit, lty = "longdash", col = col.param)
        }
      } else{
        stop(paste("the type of plot must be either 'l' (line) or 's' (stair)"))
      }
    } else if (type.cdf == "Ordinal") {
      x <- rep(cdfdata[, 1], each = 2)[-1]
      y <- rep(cdfdata[, 2], each = 2)
      tmp <- cbind(matrix(c(x, x[length(x)]), ncol = 2, byrow = TRUE), 
                   rep(NA, nrow(cdfdata)))
      x <- as.vector(t(tmp))
      tmp <- cbind(matrix(y, ncol = 2, byrow = TRUE), rep(NA, nrow(cdfdata)))
      y <- as.vector(t(tmp))
      graphics::plot(x, y, type = "l", ylim = ylimit, xlab = xlbl, ylab = ylbl,
                     ...)
      len <- length(cdfdata[, 1][tvalue])
      if (len > 1) {
        value <- rep(cdfdata[, 1][tvalue], each = 2)[-1]
        tmp <- cbind(matrix(c(value, value[length(value)]), 
                            ncol = 2, byrow = TRUE), rep(NA, len))
        value <- as.vector(t(tmp))
        len <- length(cdfdata[, 4][tvalue])
        if (len > 1) {
          lower <- rep(cdfdata[, 3][tvalue], each = 2)
          tmp <- cbind(matrix(lower, ncol = 2, byrow = TRUE), rep(NA, len))
          lower <- as.vector(t(tmp))
          upper <- rep(cdfdata[, 4][tvalue], each = 2)
          tmp <- cbind(matrix(upper, ncol = 2, byrow = TRUE), rep(NA, len))
          upper <- as.vector(t(tmp))
          graphics::lines(value, lower, lty = 3, lwd = 1.5)
          graphics::lines(value, upper, lty = 3, lwd = 1.5)
        }
      }
    } else {
      stop(paste("the type of CDF must be either 'continuous' or 'ordinal'"))
    }
    graphics::title(figlab, line = 1)
    rx <- range(graphics::par("usr")[1:2], cdfdata[, 1])
    ry <- range(graphics::par("usr")[3:4], cdfdata[, 2])
    if (legloc == "BR") {
      xjust <- 1
      yjust <- 0
      legx <- rx[2]
      legy <- ry[1]
    } else if (legloc == "BL") {
      xjust <- 0
      yjust <- 0
      legx <- rx[1]
      legy <- ry[1]
    } else if (legloc == "TR") {
      xjust <- 1
      yjust <- 1
      legx <- rx[2]
      legy <- ry[2]
    } else if (legloc == "TL") {
      xjust <- 0
      yjust <- 1
      legx <- rx[1]
      legy <- ry[2]
    }
    if (show.param) {
      mea <- round(mean(cdfdata[, 1]), round)
      mea_sd <- round(stats::sd(cdfdata[, 1]), round)
      med <- round(stats::median(cdfdata[, 1]), round)
      per <- round(spsurvey::interp.cdf(conflev, cdfest[, 6], cdfdata[, 1]), round)
      graphics::legend(x = legx, y = legy, xjust = xjust, yjust = yjust,
             legend = c("CDF Estimate",  paste(conflev, "% CL", sep = ""),
                        paste("Mean = ", mea, " (", mea_sd, ")", sep = ""),
                        paste("Median = ", med, sep = ""),
                        paste("P", conflev, " = ", per, sep = "")),
             lty = c("solid", "dotted", "dashed", "dotdash", "longdash"),
             lwd = c(1, 1.5, 1, 1, 1), bty = "n", cex = 1,
             col = c("black", "black", col.param, col.param, col.param))
    }
    else {
      graphics::legend(x = legx, y = legy, xjust = xjust, yjust = yjust,
             legend = c("CDF Estimate", paste(conflev, "% CL", sep = "")),
             lty = c(1, 3), lwd = c(1, 1.5), bty = "n", cex = 1) 
    }
    if (!is.null(ylbl.r)) {
      yl.lab <- seq(graphics::par("yaxp")[1], graphics::par("yaxp")[2], 
                    len = graphics::par("yaxp")[3] + 1)
      if (ylbl.r == "Same") {
        graphics::axis(side = 4, at = yl.lab, labels = yl.lab)
        graphics::mtext(ylbl, side = 4, line = 2, cex = graphics::par("cex"))
      } else {
        yr.lab <- spsurvey::interp.axis(yl.lab, cdfdata[, 2], cdfdata[, 5])
        graphics::axis(side = 4, at = yl.lab, 
                       labels = as.character(round(yr.lab)))
        graphics::mtext(ylbl.r, side = 4, line = 2, cex = graphics::par("cex"))
      }
    }
    graphics::par(op)
    invisible(NULL)
  }
