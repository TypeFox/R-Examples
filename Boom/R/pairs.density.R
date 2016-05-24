PairsDensity <- function(draws,
                         nlevels = 20,
                         lty = NULL,
                         color = NULL,
                         subset = NULL,
                         labels,
                         legend.location = "top",
                         legend.cex = 1,
                         label.cex = 1,
                         ...) {
  ## Produces a pairs plot showing the posterior distribution of the
  ## given list of Monte Carlo draws.  Plots above the diagonal show
  ## the posterior distribution on a scale just wide enough to fit the
  ## plots.  The diagonal shows a marginal density plot, and the sub
  ## diagonal shows the distribution with all plots on a common scale.
  ##
  ## Args:
  ##   draws: Either a matrix or a list of matrices.  If a list is
  ##     provided then each list element is plotted as a separate set
  ##     of contours, and all matrices must have the same number of
  ##     columns (though rows can differ).
  ##   nlevels:  The number of contour levels to plot.
  ##   lty:  The line types to use for the different elements in 'draws'.
  ##   color:  The color to use for different elemetns in 'draws'.
  ##   subset: If draws is a list, then this can be either a numerical
  ##     vector.  If draws has names, then subset can be a character vector.
  ##   labels: If labels is missing and the first element of draws has
  ##     non-NULL 'colnames' then these will be used to label the
  ##     pairs plot.  If a character vector of length ncol(draws[[1]])
  ##     then this character vector will be used in place of the
  ##     colnames.  If NULL then no labels will be used.
  ##   legend.location: Either "top", or "right" specifying the
  ##     location for the legend, or NULL, indicating that no legend
  ##     is desired.  if draws is a matrix or a singleton list then no
  ##     legend is produced.
  ##   legend.cex:  Scale factor to use for the legend labels.
  ##   ...: Extra arguments (graphical parameters), passed to plot,
  ##     PlotDensityContours, axis, and AddExternalLegend.

  is.odd <- function(j) {
    (j %% 2) == 1
  }

  is.even <- function(j) {
    (j %% 2) == 0
  }

  if (is.data.frame(draws)) {
    draws <- list(as.matrix(draws))
  }

  if (!is.list(draws)) {
    draws <- list(draws)
  }

  if (is.null(subset)) {
    subset <- 1:length(draws)
  }

  stopifnot(length(subset) >= 1)
  scale <- range(unlist(draws), na.rm = TRUE)

  dimension <- ncol(draws[[1]])
  number.of.groups <- length(draws)

  if (missing(labels)) {
    labels <- colnames(draws[[1]])
  }
  if (!is.null(labels)) {
    stopifnot(is.character(labels) && length(labels) == dimension)
    use.labels <- TRUE
  } else {
    ## Labels will be NULL if draws[[1]] has no colnames, or if it was
    ## supplied as NULL by the caller.  In either case is.null(labels)
    ## is the signal that determines whether labels are used.
    use.labels <- FALSE
  }

  if (number.of.groups > 1 && !is.null(legend.location)) {
    if (!is.null(names(draws))) {
      legend.labels <- names(draws)
    } else {
      legend.labels <- as.character(1:number.of.groups)
    }
    ## The following call has x.axis = FALSE because there will be no
    ## horizontal axis on the top of the plot, even though there will
    ## be one on the bottom.
    opar <- ExternalLegendLayout(dimension, dimension, legend.labels,
                                 outer.margin.lines = 4,
                                 legend.location = legend.location,
                                 legend.cex = legend.cex,
                                 x.axis = FALSE)
  } else {
    layout.matrix <- matrix(1:dimension^2, byrow = TRUE, nrow = dimension)
    opar <- par(oma = rep(4, 4), mar = rep(0, 4))
    layout(layout.matrix)
  }
  on.exit({par(opar); layout(1)})

  if (is.null(color)) {
    color <- 1:number.of.groups
  }

  if (is.null(lty)) {
    lty <- 1:number.of.groups
  }

  for (i in 1:dimension) {
    for (j in 1:dimension) {
      if (i == j) {
        ## Diagonal plot
        if (number.of.groups == 1) {
          plot(density(draws[[1]][, i]),
               axes = F,
               xlab = "",
               ylab = "",
               main = "",
               sub = "",
               ...)
        } else {
          marginal.draws <- lapply(draws, function(d) d[, i])

          CompareDensities(marginal.draws, lty = 1:number.of.groups,
                           col = color, axes = FALSE, legend.location = NULL)
        }
        box()
      } else if (j > i) {
        ## Plot is above the diagonal
        PlotDensityContours(draws,
                            j,  ## Reverse the order to go from
                            i,  ## matrix to function notation.
                            subset = subset,
                            color = color,
                            lty = lty,
                            nlevels = nlevels,
                            axes = FALSE,
                            ...)
        box()
      } else {
        ## Plot is below the diagonal.
        PlotDensityContours(draws,
                            j,
                            i,
                            subset = subset,
                            color = color,
                            lty = lty,
                            xlim = scale,
                            ylim = scale,
                            nlevels = nlevels,
                            axes = FALSE,
                            ...)
        box()
        if (i == dimension && is.odd(j)) axis(1, xpd = NA, ...)
        if (j == 1 && is.even(i)) axis(2, xpd = NA, ...)
      }
      if (use.labels) {
        if (i == 1) {
          mtext(labels[j], side = 3, cex = label.cex)
        }
        if (j == dimension) {
          mtext(labels[i], side = 4, cex = label.cex)
        }
      }
    }
  }

  if (number.of.groups > 1 && !is.null(legend.location)) {
    AddExternalLegend(legend.labels,
                      legend.location,
                      lty = lty,
                      col = color,
                      legend.cex = legend.cex,
                      ...)
  }
  return(invisible(NULL))
}

PlotDensityContours <- function(draws,
                                x.index = 1,
                                y.index = 2,
                                xlim = NULL,
                                ylim = NULL,
                                nlevels = 20,
                                subset = NULL,
                                color = NULL,
                                lty = NULL,
                                axes = TRUE,
                                ...) {
  ## Creates a contour plot of the bivariate densities in draws.  If
  ## draws is a list of Monte Carlo matrices, then each element in the
  ## list gets its own set of density contours.
  ##
  ## This is intended primarily as a panel function to use with
  ## PairsDensity.
  ##
  ## Args:
  ##   draws: Either a matrix of draws, or a list of matrices.  Each
  ##     row represents a different draw, and each column a different
  ##     variable.  If a list of matrices, each matrix should have the
  ##     same number of columns.
  ##   x.index: The index of the parameter to plot on the
  ##     horizonal axis, starting from 1.
  ##   y.index:  The index of the beta coefficient to plot on the
  ##     vertical axis, starting from 1.
  ##   xlim, ylim: If NULL then the plot will be scaled locally.
  ##     Otherwise, the limits of the horizontal and vertical axes.
  ##   nlevels:  The number of levels to use for the contour plot.
  ##   subset: The subset of the elements in draws to plot.  This can
  ##     be a numeric vector.  If draws is a named list it can also be
  ##     names.  If NULL then everything is plotted.
  ##   color: A vector of colors to use for the different subsets
  ##     begin plotted.
  ##   lty: A vector giving the type of line to use for each subset of
  ##     draws.
  ##   axes:  logical.  Should axes and a box be drawn around the figure?
  ##   ...:  Extra arguments passed to 'contour'

  if (!is.list(draws)) {
    draws <- list(draws)
  }

  stopifnot(is.list(draws))
  if (is.null(subset)) {
    subset <- 1:length(draws)
  }

  first.x <- draws[[subset[1]]][, x.index]
  first.y <- draws[[subset[1]]][, y.index]

  fixed.xlim <- !is.null(xlim)
  fixed.ylim <- !is.null(ylim)

  for (s in subset) {
    if (!fixed.xlim) {
      xlim <- range(xlim, draws[[s]][, x.index])
    }
    if (!fixed.ylim) {
      ylim <- range(ylim, draws[[s]][, y.index])
    }
  }

  if (is.null(color)) {
    color <- 1:length(subset)
  }

  if (is.null(lty)) {
    lty = 1:length(subset)
  }

  counter <- 0
  for (s in subset) {
    draws.matrix <- draws[[s]][, c(x.index, y.index)]
    counter <- counter + 1
    contour(MASS::kde2d(draws.matrix[, 1], draws.matrix[, 2]),
            nlevels = nlevels,
            xlim = xlim,
            ylim = ylim,
            col = color[counter],
            lty = lty[counter],
            add = (counter > 1),
            axes = axes,
            ...)
  }
}
