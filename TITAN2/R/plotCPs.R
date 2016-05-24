#' Plots probability densities of empirical distributions of
#' bootstrapped change points
#'
#' This function allows more detailed exploration of taxon-specific
#' response documented by TITAN through analysis of empirical
#' distributions of bootstrapped change points, comparison of those
#' distributions with observed counts, and aggregate (optionally
#' weighted) summaries of those distributions across taxa.
#'
#' Following the intitial (v1.0) TITAN publications, it was clear
#' that substantial information regarding taxon-specific change
#' points was lost when bootstrapped distributions were represented
#' solely as quantiles (i.e., as in 'plotTaxa' and the 'sppmax'
#' output table).  Empirical probability densities allow greater
#' detail and more nuanced interpretation associated with this
#' uncertainty, especially when compared against observed abundance
#' and occurrence.  Further, comparison of the summed probability
#' densities and the filtered sum(z) plots with the default sum(z)
#' output provides additional support of community changes
#' consistent with threshold behavior.
#'
#' @param titan.out A TITAN output object.
#' @param taxa.dist A logical specifying whether taxon-specific
#'   distributions should be plotted.
#' @param z.weights A logical specifying whether taxon-specific or
#'   aggregate community distributions should be weighted by their
#'   median z scores (median of z-score maxima values across
#'   bootstrap replicates).
#' @param taxaID An index specifying whether a particular taxon
#'   should be targeted for plotting. A 'NULL' value indicates all
#'   taxa should be plotted. Values >0 will select pure and
#'   reliabile taxa by their row number within the 'sppmax' output
#'   table. Character strings may also be used corresponding to the
#'   row name within the 'sppmax' output table.
#' @param cp.med A logical specifying whether change point locations
#'   should be plotted using the median value across all bootstrap
#'   replicates instead of the observed value.
#' @param cp.trace A logical specifying whether IndVals and z scores
#'   across all candidate change points should be plotted.
#' @param cp.hist A logical specifying whether histograms of
#'   replicate change point PDFs should be plotted.
#' @param stacked A logical specifying whether community level
#'   aggregations of change point PDFs are stacked or plotted
#'   separately.
#' @param xlabel A character string for the x axis label.
#' @param xmin A graphical argument specifying the value of the x
#'   axis minimum.
#' @param xmax A graphical argument specifying the value of the x
#'   axis maximum.
#' @param tck A graphical argument specifying the scale of axis tick
#'   marks.
#' @param bty A graphical argument.
#' @param ntick A graphical argument specifying the default number
#'   of axis tick marks.
#' @param cex A graphical argument specifying the scaling of the
#'   figure.
#' @param cex.axis A graphical argument specifying the scaling of
#'   the axes.
#' @param cex.leg A graphical argument specifying the scaling of the
#'   legend.
#' @param cex.lab A graphical argument specifying the scaling of the
#'   lables.
#' @param write A logical specifying whether taxa subsets are
#'   written to screen.
#' @param leg.x A graphical argument specifying the x coordinate of
#'   the legend.
#' @param leg.y A graphical argument specifying the y coordinate of
#'   the legend.
#' @param leg A logical specifying whether or not to plot the
#'   legend.
#' @param ... An argument for passing generic plotting function
#'   parameters.
#' @return Three types of plots are possible outcomes of this
#'   function.  The first (taxa.dist=T, taxID=NULL) is a matrix of
#'   histograms showing empirical distributions of bootstrapped
#'   change-point locations (as probability densities) for all pure
#'   and reliable taxa.  The value of the probability densities can
#'   be weighted by the median z score for each taxon (z.weights=T).
#'   The second plot (taxa.dist=T, taxID>0 or a taxon label)
#'   overlays a taxon-specific histogram on an abundance scatter
#'   plot and the observed change-point location.  The third plot
#'   (taxa.dist=F) shows the sum of probability densities across all
#'   pure and reliable taxa, optionally weighted by median z scores
#'   (z.weights=T) or stacked (stacked=T).
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
#' @seealso \code{\link{plotTaxa}}, \code{\link{plotSumz}}
#' @keywords TITAN kwd2
#' @export
#' @examples
#'
#' data(glades.titan)
#' plotCPs(glades.titan,
#'   taxa.dist = FALSE,
#'   xlabel = "Surface Water TP (ug/l)",
#'   stacked = TRUE
#' )
#'
#'
plotCPs <- function(titan.out, taxa.dist = T, z.weights = T, taxaID = NULL,
  cp.med = F, cp.trace = F, cp.hist = T, stacked = F, xlabel = "Environmental Gradient",
  xmin = min(titan.out$env), xmax = max(titan.out$envcls) * 1.25,
  tck = 0.025, bty = "u", ntick = 6, cex = 1.75, cex.axis = 1.75,
  cex.leg = 1.5, cex.lab = 1.75, write = F, leg.x = 0.8, leg.y = 0.8,
  leg = TRUE, ...) {

  ## PRELIMINARY CRITERIA
  minSplt <- titan.out$arguments[[1]]
  boot <- titan.out$arguments[[3]] > 0.5
  imax <- titan.out$arguments[[5]] > 0.5
  if (boot == F) {
    stop("Bootstrap Output Required for this TITAN plot")
  }
  numUnit <- length(titan.out$env)
  selTaxa = which(titan.out$sppmax[, 16] > 0)
  subsel1 = titan.out$sppmax[selTaxa, 16] == 1
  subsel2 = titan.out$sppmax[selTaxa, 16] == 2
  numTxa = length(selTaxa)


  # Blank Matrices
  cp.eden = matrix(0, numTxa, numUnit)
  cp.edenz = matrix(0, numTxa, numUnit)
  cp.eden.zmd = matrix(0, numTxa, numUnit)
  envseq = seq(from = 0, to = ceiling(max(titan.out$srtEnv)),
    by = 0.5)
  zmed = titan.out$sppmax[selTaxa, 15]

  # For each significant indicator taxa, develop weighted (by med
  # z score) and unweighted EDF of bootstrapped CPs
  for (i in 1:numTxa) {
    fe <- round(sort(titan.out$metricArray[selTaxa[i], 2, ]),
      digits = 2)
    # has quite a few non-unique values
    dens <- table(fe)/length(fe)
    edfe = approxfun(unique(fe), dens)
    for (j in 1:numUnit) {
      if (is.na(edfe(titan.out$srtEnv[j]))) {
        cp.eden[i, j] = 0
        cp.edenz[i, j] = 0
      } else {
        cp.eden[i, j] <- edfe(titan.out$srtEnv[j])
        cp.eden.zmd[i, j] <- edfe(titan.out$srtEnv[j]) *
          zmed[i]
      }
    }
  }

  ## If taxa.dist==T, plot a matrix of histograms showing
  ## frequency distributions of change points
  if (taxa.dist) {
    maxZmed = apply(cp.eden.zmd, 1, which.max)
    txaDens = matrix(NA, nrow(cp.eden), 5)
    colnames(txaDens) = c("env.max.cp", "which.env", "EDF.max",
      "z.median", "zEDF.max")
    rownames(txaDens) = rownames(titan.out$sppmax[selTaxa,
      ])
    for (i in 1:nrow(cp.eden)) {
      txaDens[i, 1] = round((titan.out$srtEnv[maxZmed[i] -
        1] + titan.out$srtEnv[maxZmed[i]])/2, digits = 2)
      txaDens[i, 2] = maxZmed[i]
      txaDens[i, 3] = round(cp.eden[i, maxZmed[i]], digits = 2)
      txaDens[i, 4] = round(zmed[i], digits = 2)
      txaDens[i, 5] = round(cp.eden.zmd[i, maxZmed[i]], digits = 2)
    }

    if (!(is.null(taxaID))) {
      if (is.numeric(taxaID)) {
        taxNum = taxaID
      } else {
        taxNum = which(taxaID == rownames(titan.out$sppmax))
      }
      if (titan.out$sppmax[taxNum, 16] < 1) {
        stop("This taxon is either impure or unreliable, try a robust taxon")
      }
      tag = rep(0, nrow(titan.out$sppmax))
      tag[taxNum] <- 1
      tagNum <- which(tag[selTaxa] > 0)
      par(mar = c(5, 5, 4, 2), oma = c(0, 0, 0, 3))
      plot(titan.out$env, titan.out$taxa[, taxNum], xlab = xlabel,
        ylab = "Abundance", axes = T, col = "black", cex.axis = cex.axis,
        cex = cex, cex.lab = cex.lab, tck = tck)
      segments(titan.out$sppmax[taxNum, 8], max(titan.out$taxa[,
        taxNum], na.rm = T), titan.out$sppmax[taxNum, 12],
        max(titan.out$taxa[, taxNum], na.rm = T), col = "red",
        lwd = 2)
      if (cp.med) {
        cp.choice = titan.out$sppmax[taxNum, 10]
      } else {
        if (imax) {
          cp.choice = titan.out$sppmax[taxNum, 1]
        } else {
          cp.choice = titan.out$sppmax[taxNum, 2]
        }
      }
      symbols(cp.choice, max(titan.out$taxa[, taxNum], na.rm = T),
        circles = titan.out$sppmax[taxNum, 15], inches = 0.1,
        add = TRUE, fg = "red", bg = "white", lwd = 2)
      if (cp.trace) {
        par(new = T)
        plot(titan.out$srtEnv, c(rep(NA, minSplt), titan.out$ivzScores[nrow(titan.out$sppmax) +
          taxNum, ], rep(NA, minSplt - 1)), type = "l",
          axes = F, xlab = "", ylab = "", col = "red",
          lty = 2, ylim = c(0, max(titan.out$ivzScores[nrow(titan.out$sppmax) +
          taxNum, ], na.rm = T)))

        par(new = T)
        plot(titan.out$srtEnv, c(rep(NA, minSplt), titan.out$ivzScores[(nrow(titan.out$sppmax) *
          2) + taxNum, ], rep(NA, minSplt - 1)), type = "l",
          axes = F, xlab = "", ylab = "", col = "grey80",
          lty = 2, ylim = c(0, max(titan.out$ivzScores[(nrow(titan.out$sppmax) *
          2) + taxNum, ], na.rm = T)))
      }
      if (cp.hist) {
        par(new = T)
        plot(titan.out$srtEnv, cp.eden[tagNum, ], type = "h",
          axes = F, xlab = "", ylab = "", col = "blue",
          main = rownames(titan.out$sppmax)[taxNum], ylim = c(0,
          max(cp.eden[tagNum, ], na.rm = T)))
        axis(4, pretty(c(0, max(cp.eden[tagNum, ], na.rm = T)),
          6), cex.axis = cex.axis, tck = tck, mgp = c(2.5,
          0.5, 0))
        mtext("Density", side = 4, line = 3, cex = cex.lab)
      }
    } else {

      par(mfrow = c(ceiling(sqrt(length(selTaxa))), ceiling(sqrt(length(selTaxa)))),
        mar = c(2, 2, 2, 2))
      ## If z.weights==T, densities are weighted by median z scores
      if (z.weights) {
        if (sum(subsel1, na.rm = T) > 0) {
          for (i in 1:nrow(cp.eden[subsel1, ])) {
          plot(titan.out$srtEnv, cp.eden.zmd[subsel1,
            ][i, ], type = "h", axes = T, xlab = "",
            ylab = "", col = "blue", main = rownames(titan.out$sppmax)[selTaxa[subsel1]][i],
            ylim = c(0, 1.5))
          }
        }
        if (sum(subsel2, na.rm = T) > 0) {
          for (i in 1:nrow(cp.eden[subsel2, ])) {
          plot(titan.out$srtEnv, cp.eden.zmd[subsel2,
            ][i, ], type = "h", axes = T, xlab = "",
            ylab = "", col = "red", main = rownames(titan.out$sppmax)[selTaxa[subsel2]][i],
            ylim = c(0, 1.5))
          }
        }
        ## If z.weights==F, densities are raw frequencies
      } else {
        if (sum(subsel1, na.rm = T) > 0) {
          for (i in 1:nrow(cp.eden[subsel1, ])) {
          plot(titan.out$srtEnv, cp.eden[subsel1, ][i,
            ], type = "h", axes = T, xlab = "", ylab = "",
            col = "blue", main = rownames(titan.out$sppmax)[selTaxa[subsel1]][i],
            ylim = c(0, 0.3))
          }
        }
        if (sum(subsel2, na.rm = T) > 0) {
          for (i in 1:nrow(cp.eden[subsel2, ])) {
          plot(titan.out$srtEnv, cp.eden[subsel2, ][i,
            ], type = "h", axes = T, xlab = "", ylab = "",
            col = "red", main = rownames(titan.out$sppmax)[selTaxa[subsel2]][i],
            ylim = c(0, 0.3))
          }
        }
      }

      if (sum(subsel1, na.rm = T) > 0 & write) {
        print("Group 1 Pure and Reliable Taxa")
        print(txaDens[subsel1, ])
      }
      if (sum(subsel2, na.rm = T) > 0 & write) {
        print("Group 2 Pure and Reliable Taxa")
        print(txaDens[subsel2, ])
      }
    }

    ## If taxa.dist==F, sum the densities across all taxa
  } else {

    sumDen1 = colSums(cp.eden[subsel1, ], na.rm = T)
    sumDen2 = colSums(cp.eden[subsel2, ], na.rm = T)
    sumDen = colSums(cp.eden, na.rm = T)
    sumDenzd1 = colSums(cp.eden.zmd[subsel1, ], na.rm = T)
    sumDenzd2 = colSums(cp.eden.zmd[subsel2, ], na.rm = T)
    sumDenzd = colSums(cp.eden.zmd, na.rm = T)

    sumTab = matrix(NA, 3, 4)
    rownames(sumTab) = c("Group 1", "Group 2", "P&R Taxa")
    colnames(sumTab) = c("uw.env", "uw.max", "zw.env", "zw.max")
    maxDen1 = which.max(sumDen1)
    maxDen2 = which.max(sumDen2)
    maxDen = which.max(sumDen)
    maxZden1 = which.max(sumDenzd1)
    maxZden2 = which.max(sumDenzd2)
    maxZden = which.max(sumDenzd)

    if (sum(subsel1, na.rm = T) > 0) {
      if (maxDen1 > 1) {
        sumTab[1, 1] = round((titan.out$srtEnv[maxDen1 -
          1] + titan.out$srtEnv[maxDen1])/2, digits = 2)
      } else {
        sumTab[1, 1] = titan.out$srtEnv[maxDen1]
      }

      sumTab[1, 2] = round(sumDen1[maxDen1], digits = 2)

      if (maxDen1 > 1) {
        sumTab[1, 3] = round((titan.out$srtEnv[maxZden1 -
          1] + titan.out$srtEnv[maxZden1])/2, digits = 2)
      } else {
        sumTab[1, 3] = titan.out$srtEnv[maxZden1]
      }

      sumTab[1, 4] = round(sumDenzd1[maxZden1], digits = 2)
    }
    if (sum(subsel2, na.rm = T) > 0) {
      sumTab[2, 1] = round((titan.out$srtEnv[maxDen2 - 1] +
        titan.out$srtEnv[maxDen2])/2, digits = 2)
      sumTab[2, 2] = round(sumDen2[maxDen2], digits = 2)
      sumTab[2, 3] = round((titan.out$srtEnv[maxZden2 - 1] +
        titan.out$srtEnv[maxZden2])/2, digits = 2)
      sumTab[2, 4] = round(sumDenzd2[maxZden2], digits = 2)
    }

    if (maxDen > 1) {
      sumTab[3, 1] = round((titan.out$srtEnv[maxDen - 1] +
        titan.out$srtEnv[maxDen])/2, digits = 2)
    } else {
      sumTab[3, 1] = titan.out$srtEnv[maxDen]
    }

    sumTab[3, 2] = round(sumDen[maxDen], digits = 2)

    if (maxZden > 1) {
      sumTab[3, 3] = round((titan.out$srtEnv[maxZden - 1] +
        titan.out$srtEnv[maxZden])/2, digits = 2)
    } else {
      sumTab[3, 3] = titan.out$srtEnv[maxZden]
    }

    sumTab[3, 4] = round(sumDenzd[maxZden], digits = 2)
    print("Summary of Summed Density Functions")
    print(sumTab)

    ## Define Stacked Plot Function
    plot.stacked <- function(x, y, order.method = "as.is",
      ylab = "", xlab = "", border = NULL, lwd = 0.5, col = c("blue",
        "red"), ylim = NULL, ...) {
      if (sum(y < 0) > 0)
        warning("Y cannot contain negative numbers")
      if (is.null(border))
        border <- par("fg")
      border <- as.vector(matrix(border, nrow = ncol(y),
        ncol = 1))
      col <- as.vector(matrix(col, nrow = ncol(y), ncol = 1))
      lwd <- as.vector(matrix(lwd, nrow = ncol(y), ncol = 1))
      if (order.method == "max") {
        ord <- order(apply(y, 2, which.max))
        y <- y[, ord]
        col <- col[ord]
        border <- border[ord]
      }
      top.old <- x * 0
      polys <- vector(mode = "list", ncol(y))
      for (i in seq(polys)) {
        top.new <- top.old + y[, i]
        polys[[i]] <- list(x = c(x, rev(x)), y = c(top.old,
          rev(top.new)))
        top.old <- top.new
      }
      if (is.null(ylim))
        ylim <- range(sapply(polys, function(x) range(x$y,
          na.rm = TRUE)), na.rm = TRUE)
      plot(x, y[, 1], ylab = ylab, xlab = xlab, ylim = ylim,
        t = "n", col = "blue", ...)
      for (i in seq(polys)) {
        polygon(polys[[i]], border = col[i], col = col[i],
          lwd = lwd[i])
      }
    }

    ## If z.weights==T, summed densities are weighted by median z
    ## scores
    par(mar = c(5, 5, 4, 2), oma = c(0, 0, 0, 3))
    if (z.weights) {
      if (stacked) {
        plot.stacked(titan.out$srtEnv, cbind(sumDenzd1,
          sumDenzd2), xlab = xlabel, ylab = "Summed z-Weighted Probability Densities",
          sub = "pure and reliable taxa only", xlim = c(xmin,
          xmax), cex = 1.75, cex.axis = 1.75, tck = 0.025,
          cex.lab = 1.5)
        points(titan.out$srtEnv[which.max(sumDenzd1)],
          max(sumDenzd1, na.rm = T), pch = 19, col = "blue")
        points(titan.out$srtEnv[which.max(sumDenzd2)],
          max(sumDenzd2, na.rm = T), pch = 19, col = "red")
      } else {
        plot(titan.out$srtEnv, sumDenzd, xlab = xlabel,
          ylab = "Summed z-Weighted Probability Densities",
          ty = "b", col = "gray", lty = 1, sub = "pure and reliable taxa only",
          xlim = c(xmin, xmax), cex = 1.75, cex.axis = 1.75,
          tck = 0.025, cex.lab = 1.75)
        points(titan.out$srtEnv, sumDenzd2, ty = "b", col = "red",
          pch = 19)
        points(titan.out$srtEnv, sumDenzd1, ty = "b", col = "blue",
          pch = 19)
        # points(titan.out$srtEnv,sumDenzd,ty='b',col='gray')
      }
      ## If z.weights==F, densities are unweighted
    } else {
      if (stacked) {
        plot.stacked(titan.out$srtEnv, cbind(sumDen1, sumDen2),
          xlab = xlabel, ylab = "Summed Unweighted Probability Densities",
          sub = "pure and reliable taxa only", xlim = c(xmin,
          xmax), cex = 1.75, cex.axis = 1.75, tck = 0.025,
          cex.lab = 1.5)
        points(titan.out$srtEnv[which.max(sumDen1)], max(sumDen1,
          na.rm = T), pch = 19, col = "blue")
        points(titan.out$srtEnv[which.max(sumDen2)], max(sumDen2,
          na.rm = T), pch = 19, col = "red")
      } else {
        plot(titan.out$srtEnv, sumDen, xlab = xlabel, ylab = "Summed Unweighted Probability Densities",
          ty = "b", col = "gray", lty = 1, sub = "pure and reliable taxa only",
          xlim = c(xmin, xmax), cex = 1.75, cex.axis = 1.75,
          tck = 0.025, cex.lab = 1.75)
        points(titan.out$srtEnv, sumDen2, ty = "b", col = "red",
          pch = 19)
        points(titan.out$srtEnv, sumDen1, ty = "b", col = "blue",
          pch = 19)
        # points(titan.out$srtEnv,sumDen,ty='b',col='gray')
      }
    }

  }

  if (leg) {
    if (stacked) {
      legend("topright", legend = c("z- max", "z+ max"),
        col = c("blue", "red"), pch = 19)
    } else {
      if (!taxa.dist) {
        legend("topright", legend = c("all z", "z-", "z+"),
          col = c("gray", "blue", "red"), pch = 19)
      }
    }
  }
  ## End plotCPs
}
