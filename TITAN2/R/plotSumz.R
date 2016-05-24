#' Plots the pattern of community change along an environmental
#' gradient
#'
#' Creates a plot of community-level sums of taxon-specific change
#' along an environmental gradient and optionally conveys
#' uncertainty associated with the maximum community change derived
#' from decreasers or increasers.
#'
#' This function assumes that the TITAN object contains bootstrap
#' summaries and filtering information and automatically determines
#' whether this is the case.  Without the bootstrap summaries, only
#' unfiltered change magnitudes are plotted.
#'
#' The original sum(z) plots (v1.0) did not filter taxa using purity
#' and reliability.  Because these taxa often have small z scores,
#' they are unlikely to contribute significantly to the sum(z)
#' profiles.  However, subsequent investigation has demonstrated
#' that when sufficient numbers of taxa are involved, it is possible
#' for noisy data to generate artifactual peaks in low-magnitude
#' sum(z) profiles or plateaus.  Therefore, we recommend evaluating
#' filtered versions of the sum(z) to assess this potential in v2.0.
#'
#' @param titan.out A TITAN output object.
#' @param filter A logical indicating whether only pure and reliable
#'   taxa should be used to create the plot.This is the recommended
#'   as a check of the unfiltered default to assess whether impure
#'   or unreliable taxa are substantially contributing to the
#'   distribution of sum(z) scores.
#' @param cumfrq A logical specifying whether cumulative frequencies
#'   of sum(z) maxima across bootstrap replicates should be plotted.
#' @param bootz1 A logical specifying whether decreasing cumulative
#'   frequencies exist or should be plotted.
#' @param bootz2 A logical specifying whether increasing cumulative
#'   frequencies exist or should be plotted.
#' @param sumz1 A logical specifying whether decreasing changes
#'   should be plotted.
#' @param sumz2 A logical specifying whether increasing changes
#'   should be plotted.
#' @param xmin A graphical argument specifying the value of the x
#'   axis minimum.
#' @param xmax A graphical argument specifying the value of the x
#'   axis maximum.
#' @param xlabel A character string for the x axis label.
#' @param y1label A character specifying the label of the second y
#'   axis
#' @param y2label A character specifying the label of the second y
#'   axis
#' @param log A graphical argument specifying whether an axis should
#'   be log scaled.
#' @param at A graphical argument for controling placement of the x
#'   axis label
#' @param tck A graphical argument specifying the scale of axis tick
#'   marks.
#' @param bty A graphical argument.
#' @param ntick A graphical argument specifying the default number
#'   of axis tick marks.
#' @param prtty A logical specifying whether pretty() should be used
#'   to plot axis labels.
#' @param dig A numeric argument specifying the number of decimal
#'   digits in axes.
#' @param cex A graphical argument specifying the scaling of the
#'   figure.
#' @param cex.axis A graphical argument specifying the scaling of
#'   the axes.
#' @param cex.leg A graphical argument specifying the scaling of the
#'   legend.
#' @param cex.lab A graphical argument specifying the scaling of the
#'   lables.
#' @param leg.x A graphical argument specifying the x coordinate of
#'   the legend.
#' @param leg.y A graphical argument specifying the y coordinate of
#'   the legend.
#' @param legend A logical specifying whether or not to plot the
#'   legend.
#' @param pch1 A graphical argument specifying the type of group 1
#'   symbols.
#' @param pch2 A graphical argument specifying the type of group 2
#'   symbols.
#' @param col1 A graphical argument specifying the color of group 1
#'   symbols.
#' @param col2 A graphical argument specifying the color of group 2
#'   symbols.
#' @param ... An argument for passing generic plotting function
#'   parameters.
#' @return A plot of sum(z-) and sum(z+) profiles along the
#'   environmental gradient.
#' @references Baker, ME and RS King.  2010. A new method for
#'   detecting and interpreting biodiversity and ecological
#'   community thresholds.  Methods in Ecology and Evolution 1(1):
#'   25:37.
#' @references King, RS and ME Baker  2010. Considerations for
#'   identifying and interpreting ecological community thresholds.
#'   Journal of the North American Benthological Association
#'   29(3):998-1008.
#' @note Should not be used with output objects from TITAN v1.0.
#' @author M. Baker and R. King
#' @seealso \code{\link{plotTaxa}}, \code{\link{plotCPs}}
#' @keywords TITAN sum(z)
#' @export
#' @examples
#'
#' data(glades.titan)
#'
#'
plotSumz <- function(titan.out, filter = F, cumfrq = T, bootz1 = T,
  bootz2 = T, sumz1 = T, sumz2 = T, xmin = min(titan.out$env),
  xmax = max(titan.out$envcls) * 1.25, xlabel = "Environmental Gradient",
  y1label = NULL, y2label = "Cumulative Frequency", log = "",
  at = NULL, tck = 0.025, bty = "u", ntick = 6, prtty = T, dig = 1,
  cex = 1.75, cex.axis = 1.75, cex.leg = 1.5, cex.lab = 1.75,
  leg.x = 0.8, leg.y = 0.8, legend = TRUE, pch1 = 16, pch2 = 1,
  col1 = "black", col2 = "black", ...) {


  ## SET GRAPH DIMENSIONS
  par(mar = c(8, 5, 4, 2), oma = c(0, 0, 0, 3))

  if (filter) {
    ivz <- titan.out$ivz.f
    maxPsumz <- titan.out$maxFsumz
    if (is.null(y1label))
      y1label <- "Filtered Sum(z)"
  } else {
    ivz <- titan.out$ivz
    maxPsumz <- titan.out$maxSumz
    if (is.null(y1label))
      y1label <- "Unfiltered Sum(z)"
  }
  boot <- titan.out$arguments[[3]] > 0.5


  if (sum(maxPsumz[, 1], na.rm = T) < 1) {
    sumz1 = F
  }
  if (sum(maxPsumz[, 2], na.rm = T) < 1) {
    sumz2 = F
  }

  if (boot & cumfrq) {
    ## PLOT CUMULATIVE PROBABILITY CURVE FOR GROUP 1
    if (sumz1) {
      freq1 <- table(maxPsumz[, 1])
      freq2 <- table(maxPsumz[, 2])
      if (bootz1) {
        if (length(unique(maxPsumz[, 1])) > length(freq1)) {
          fil = matrix(0, length(unique(maxPsumz[, 1])) -
          length(freq1) + 1, 1)
          plot(rbind(min(maxPsumz[, 1]), matrix(sort(unique(maxPsumz[,
          1])))), rbind(fil, matrix(cumsum(freq1)/sum(freq1))),
          type = "l", lty = 1, lwd = 2, col = col1, log = log,
          axes = FALSE, ylim = c(0, 1), xlim = c(xmin,
            xmax), cex = cex.axis, xlab = "", ylab = "")
        } else {
          plot(rbind(min(maxPsumz[, 1]), matrix(sort(unique(maxPsumz[,
          1])))), rbind(0, matrix(cumsum(freq1)/sum(freq1))),
          type = "l", lty = 1, lwd = 2, col = col1, log = log,
          axes = FALSE, ylim = c(0, 1), xlim = c(xmin,
            xmax), cex = cex.axis, xlab = "", ylab = "")
        }
      }
    }

    ## PLOT CUMULATIVE PROBABILITY CURVE FOR GROUP 2
    if (sumz2) {
      if (bootz1) {
        if (bootz2) {
          if (length(unique(maxPsumz[, 2])) > length(freq2)) {
          fil2 = matrix(0, length(unique(maxPsumz[, 2])) -
            length(freq2) + 1, 1)
          points(rbind(min(maxPsumz[, 2]), matrix(sort(unique(maxPsumz[,
            2])))), rbind(fil2, matrix(cumsum(freq2)/sum(freq2))),
            type = "l", lty = 2, lwd = 2, col = col2)
          } else {
          points(rbind(min(maxPsumz[, 2]), matrix(sort(unique(maxPsumz[,
            2])))), rbind(0, matrix(cumsum(freq2)/sum(freq2))),
            type = "l", lty = 2, lwd = 2, col = col2)
          }
        }
      } else {
        if (bootz2) {
          if (length(unique(maxPsumz[, 2])) > length(freq2)) {
          fil2 = matrix(0, length(unique(maxPsumz[, 2])) -
            length(freq2) + 1, 1)
          plot(rbind(min(maxPsumz[, 2]), matrix(sort(unique(maxPsumz[,
            2])))), rbind(fil2, matrix(cumsum(freq2)/sum(freq2))),
            type = "l", lty = 2, lwd = 2, col = col2,
            log = log, axes = FALSE, ylim = c(0, 1),
            xlim = c(xmin, xmax), cex = cex.axis, xlab = "",
            ylab = "")
          } else {
          plot(rbind(min(maxPsumz[, 2]), matrix(sort(unique(maxPsumz[,
            2])))), rbind(0, matrix(cumsum(freq2)/sum(freq2))),
            type = "l", lty = 2, lwd = 2, col = col2,
            log = log, axes = FALSE, ylim = c(0, 1),
            xlim = c(xmin, xmax), cex = cex.axis, xlab = "",
            ylab = "")
          }
        }
      }
    }

    ## ADD SECOND Y-AXIS WITH PROBABILITY LABELS
    if (bootz1 | bootz2) {
      axis(4, pretty(c(0, 1), 6), las = 1, cex.axis = cex.axis,
        tck = tck, mgp = c(2.5, 0.5, 0))
      mtext(y2label, side = 4, line = 3, cex = cex.lab)
    }

    ## FIND MINIMUM AND MAX VALUES TO COMPLETE CURVES

    sumz1.min <- min(maxPsumz[, 1])
    sumz2.min <- min(maxPsumz[, 2])
    sumz1.max <- max(maxPsumz[, 1])
    sumz2.max <- max(maxPsumz[, 2])


    ## COMPLETE PROBABILITY CURVES
    if (bootz1) {
      segments(sumz1.max, 1, xmax, 1, col = col1, lty = 1,
        lwd = 2)
      segments(xmin, 0, sumz1.min, 0, col = col1, lty = 1,
        lwd = 2)
    }

    if (bootz2) {
      segments(xmin, 0, sumz2.min, 0, col = col2, lty = 2,
        lwd = 2)
      segments(sumz2.max, 1, xmax, 1, col = col2, lty = 2,
        lwd = 2)
    }
  }

  ## SET PLOT WINDOW TO ACCEPT NEW PARAMETERS
  if (boot & cumfrq) {
    if (bootz1 | bootz2) {
      par(new = TRUE)
    }
  }

  ## PLOT OBSERVED SUMZ SCORES FOR BOTH GROUPS AND ADD APPROPRIATE
  ## LEGEND
  if (sumz1) {
    plot(titan.out$envcls, ivz[, 1], type = "b", pch = pch1,
      col = col1, lty = 1, xlim = c(xmin, xmax), ylim = c(min(ivz,
        na.rm = T), max(ivz, na.rm = T)), log = log, xlab = xlabel,
      cex.lab = cex.lab, axes = FALSE, ylab = "")
    if (sumz2) {
      points(titan.out$envcls, ivz[, 2], type = "b", col = col2,
        lty = 2, pch = pch2)
    }
  } else {
    if (sumz2) {
      plot(titan.out$envcls, ivz[, 2], type = "b", col = col2,
        lty = 1, pch = pch2, xlim = c(xmin, xmax), ylim = c(min(ivz,
          na.rm = T), max(ivz, na.rm = T)), log = log,
        xlab = xlabel, cex.lab = cex.lab, axes = FALSE,
        ylab = "")
    }
  }

  ## ADD X-AXIS WITH SAME RANGE AS DUMMY PLOT
  if (log == "x") {
    axis(1, at = at, mgp = c(2, 0.5, 0), cex.axis = cex.axis,
      tck = tck)
  } else {
    if (prtty) {
      axis(1, pretty(xmin:xmax, ntick), mgp = c(2, 0.5, 0),
        cex.axis = cex.axis, tck = tck)
    } else {
      axis(1, at = seq(from = round(xmin, digits = dig),
        to = round(xmax, digits = dig), by = round((xmax -
          xmin)/4, digits = dig)), mgp = c(2, 0.5, 0),
        cex.axis = cex.axis, tck = tck)
    }
  }

  axis(2, pretty(min(ivz, na.rm = T):max(ivz, na.rm = T), 6),
    mgp = c(5, 0.5, 0), las = 1, cex.axis = cex.axis, tck = tck)
  mtext(y1label, side = 2, line = 3.5, cex = cex.lab)

  ## ADD BOX AND LEGEND
  if (bootz1) {
    box(which = "plot", bty = bty)
  } else {
    if (bootz2) {
      box(which = "plot", bty = bty)
    }
  }

  if (sumz1) {
    leg = c("z-")
  } else {
    if (sumz2) {
      leg = c("z+")
    }
  }
  if (sumz1) {
    if (sumz2) {
      leg = c("z-", "z+")
    }
  }

  fill.leg = c(col1, ifelse(col2 == "black", "white", col2))
  if (legend) {
    legend(titan.out$envcls[length(titan.out$envcls) * leg.x],
      max(ivz, na.rm = T) * leg.y, bty = "n", leg, fill = fill.leg,
      ncol = 1, plot = TRUE, cex = cex.leg)
  }
  ### END plotSumz
}
