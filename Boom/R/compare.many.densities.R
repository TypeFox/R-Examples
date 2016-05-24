CompareManyDensities <- function(list.of.arrays,
                                 style = c("density", "box"),
                                 main = "",
                                 color = NULL,
                                 gap = 0,
                                 burn = 0,
                                 suppress.labels = FALSE,
                                 x.same.scale = TRUE,
                                 y.same.scale = FALSE,
                                 xlim = NULL,
                                 ylim = NULL,
                                 legend.location = c("top", "right"),
                                 legend.cex = 1,
                                 reflines = NULL,
                                 ...) {
  ## Produce a plot that compares the marginal distribution of each
  ## element of a vector or matrix across several groups.
  ##
  ## Args:
  ##   list.of.arrays: A list of arrays representing the MCMC draws of
  ##     the vector or matrix in question.  Each list element
  ##     represents a different group.  The first index in each list
  ##     list element represents the Monte Carlo draw number (or
  ##     iteration).  The remaining indices represent the variables to
  ##     be plotted.  If the first list element has variable names
  ##     assigned to its indices, these will be used to label the
  ##     plots.
  ##   style:  The style of plot to use for comparing distributions.
  ##   main:  The main title of the plot.
  ##   color:  A vector of colors to be used for representing the groups.
  ##   gap:  The gap (in lines) between plots).
  ##   burn:  The number of MCMC iterations to be discarded as burn-in.
  ##   suppress.labels: Logical.  If FALSE then the dimnames (if any)
  ##     of the first element in list.of.arrays will be used to
  ##     annotate the plot.  If TRUE then no labels will be used.
  ##   x.same.scale: Logical indicating whether the same horizontal
  ##     scale should be used for all the plots.
  ##   y.same.scale: Logical indicating whether the same vertical
  ##     scale should be used for all the plots.  This argument is
  ##     ignored if style == "box".
  ##   xlim: Either NULL, or a pair of numbers giving limits for the
  ##     horizontal axis.  If xlim is set then the same xlim values
  ##     will be used for all plots and the x.same.scale argument will
  ##     be ignored.
  ##   ylim: Either NULL, or a pair of numbers giving limits for the
  ##     vertical axis.  If ylim is set then the same ylim values
  ##     will be used for all plots and the y.same.scale argument will
  ##     be ignored.  This argument is ignored if style == "box".
  ##   legend.location: The location of the legend, either on top or
  ##     at the right.  It can also be NULL in which case no legend
  ##     will appear.  The legend names will be taken from
  ##     names(list.of.arrays).  If it does not have names, then no
  ##     legend will be produced.
  ##   legend.cex:  The relative scale factor to use for the legend text.
  ##   reflines: This can be NULL, in which case no reference lines
  ##     are drawn, it can be a single real number in which case a
  ##     reference line will be drawn at that value in each panel, or
  ##     it can be a vector with length equal to the number of panels,
  ##     in which case a reference line will be drawn at each
  ##     panel-specific value.
  ##   ...: Extra arguments passed to either CompareDensities, or
  ##     boxplot.
  ##
  GetDensityLimits <- function(x) {
    ## Args:
    ##   x is a matrix or a 3-way array
    ## Returns:
    ##   The upper limit of the density() function applied to each
    ##   column in x.
    dmax <- function(x) {return(max(density(x)$y))}
    if (is.matrix(x)) {
      return(max(apply(x, 2, dmax)))
    } else {
      return(max(apply(x, c(2, 3), dmax)))
    }
  }

  style <- match.arg(style)
  stopifnot(is.list(list.of.arrays))
  stopifnot(all(sapply(list.of.arrays, function(x) length(dim(x))) ==
  length(dim(list.of.arrays[[1]]))))

  if (!is.null(reflines)) stopifnot(is.numeric(reflines))

  ## Deduce the number of rows and column necessary, and find row,
  ## column, or plot labels.
  ndim <- length(dim(list.of.arrays[[1]]) )
  have.matrix.labels <- FALSE
  vector.labels <- NULL
  if (ndim == 3) {
    ## Handle the case of matrix valued inputs.  The first index
    ## counts Monte Carlo draws.
    dims <- dim(list.of.arrays[[1]])
    nr <- dims[2]
    nc <- dims[3]
    nvars <- nr * nc
    matrix.data <- TRUE
    dim.names <- dimnames(list.of.arrays[[1]])
    if (!is.null(dim.names)) {
      have.matrix.labels <- TRUE
      row.labels <- dim.names[[2]]
      column.labels <- dim.names[[3]]
    }
  } else if (ndim == 2) {
    ## Handle the case of vector valued inputs.  The first index
    ## counts Monte Carlo draws.
    nvars <- ncol(list.of.arrays[[1]])
    nr <- floor(sqrt(nvars))
    nc <- ceiling(nvars / nr)
    matrix.data <- FALSE
    vector.labels <- colnames(list.of.arrays[[1]])
  } else {
    stop("The input to CompareManyDensities must be a list of ",
         "matrices or 3-way arrays.")
  }

  layout.matrix <- matrix(1:(nr * nc), nrow = nr, ncol = nc, byrow = TRUE)
  if (!is.null(legend.location) && !is.null(names(list.of.arrays))) {
    use.legend <- TRUE
    legend.location <- match.arg(legend.location)
    number.of.plots <- nr * nc
    if (legend.location == "top") {
      layout.matrix <- rbind(number.of.plots + 1, layout.matrix)
    } else if (legend.location == "right") {
      layout.matrix <- cbind(layout.matrix, number.of.plots + 1)
    }
  } else {
    use.legend <- FALSE
  }

  ## Set graphical parameters: margins and the number of plots.  If a
  ## title was supplied we need to leave additional space to print it.
  top.margin <- 4
  if (main != "") {
    top.margin <- top.margin + 4.1
  }
  if (use.legend) {
    if (legend.location == "top") {
      legend.buffer <- max(strheight(names(list.of.arrays),
                                     units = "inches",
                                     cex = legend.cex))
      legend.buffer <- legend.buffer * 1 * 2.54
      if (x.same.scale) {
        legend.buffer <- legend.buffer + .5
      }
      layout(layout.matrix,
             heights = c(lcm(legend.buffer), rep(1, nr)))
    } else if (legend.location == "right") {
      legend.buffer <- max(strwidth(names(list.of.arrays),
                                    units = "inches",
                                    cex = legend.cex))
      legend.buffer <- legend.buffer * 2.54 + legend.cex
      layout(layout.matrix, widths = c(rep(1, nc), lcm(legend.buffer)))
    }
  } else {
    layout(layout.matrix)
  }

  original.par <- par(oma = c(4, 4, top.margin, 4),
                      mar = rep(gap/2, 4),
                      xpd = FALSE)
  on.exit(par(original.par))

  ## Deduce xlim and ylim, if necessary.
  if (x.same.scale && is.null(xlim)) {
    xlim <- range(unlist(list.of.arrays), na.rm = TRUE)
  }
  if (style == "box") {
    ylim <- NULL
  }
  if (y.same.scale && is.null(ylim) && style == "density") {
    ylim <- range(c(0, sapply(list.of.arrays, GetDensityLimits)),
                  na.rm = TRUE)
  }

  ## Set up the color scheme, if necessary.
  number.of.groups <- length(list.of.arrays)
  if (is.null(color)) {
    color <- 1:number.of.groups + 1 ## skip black
  }

  variable.number <- 0
  for (i in 1:nr) {
    for (j in 1:nc) {
      variable.number <- variable.number + 1
      if (variable.number <= nvars) {
        ## Assemble the list of values to be plotted in this panel.
        values <- list()
        for (g in 1:number.of.groups) {
          iterations <- 1:(dim(list.of.arrays[[g]])[1])
          if (burn > 0) {
            iterations <- iterations[-(1:burn)]
          }
          if (matrix.data) {
            values[[g]] <- list.of.arrays[[g]][iterations, i, j]
          } else {
            values[[g]] <- list.of.arrays[[g]][iterations, variable.number]
          }
        }

        ## Make the plot for this panel.
        if (style == "density") {
          CompareDensities(values, col = color, xlim = xlim, ylim = ylim,
                           xaxt = "n", yaxt = "n", ...)
          if (!is.null(reflines)) {
            if (length(reflines) == 1) {
              abline(v = reflines, lwd = 2, ...)
            } else {
              abline(v = reflines[variable.number], lwd = 2, ...)
            }
          }
        } else if (style == "box") {
          boxplot(values, col = color, xaxt = "n", yaxt = "n",
                  ylim = xlim, horizontal = TRUE, ...)
          if (!is.null(reflines)) {
            if (length(reflines) == 1) {
              abline(h = reflines, lwd = 2, ...)
            } else {
              abline(h = reflines[variable.number], lwd = 2, ...)
            }
          }
        }

        ## If the panel needs a label, add the label.
        if (!suppress.labels && !is.null(vector.labels)) {
          mtext(vector.labels[variable.number], line = -2)
        }

      } else {
        ## This block puts empty plots on the bottom row that are
        ## place holders for shared axes.
        if (is.null(ylim)) {
          ylim <- range(unlist(list.of.arrays), na.rm = TRUE)
        }
        if (is.null(xlim)) {
          xlim <- range(unlist(list.of.arrays), na.rm = TRUE)
        }
        plot(xlim, ylim, type = "n", xaxt = "n", yaxt = "n")
      }

      ##------------------------------------------------------------
      ## Put axes and row/column labels where they belong.
      if (i == nr) {
        ## Bottom row
        if (IsOdd(j)) {
          if (x.same.scale) {
            axis(1)
          }
        } else {
          if (have.matrix.labels && !suppress.labels) {
            mtext(column.labels[j], side = 1, line = 2)
          }
        }
      }

      if (i == 1) {
        ## top row
        if (IsEven(j)) {
          if (x.same.scale) {
            axis(3)
          }
        } else {
          if (have.matrix.labels && !suppress.labels) {
            mtext(column.labels[j], side = 3, line = 2)
          }
        }
      }

      if (j == nc) {
        ## right column
        if (IsEven(i)) {
          if (y.same.scale && style == "density") {
            axis(4)
          }
        } else {
          if (have.matrix.labels && !suppress.labels) {
            mtext(row.labels[i], side = 4, line = 2)
          }
        }
      }

      if (j == 1) {
        ## left column
        if (IsOdd(i)) {
          if (y.same.scale && style == "density") {
            axis(2)
          }
        } else {
          if (have.matrix.labels && !suppress.labels) {
            mtext(row.labels[i], side = 2, line = 2)
          }
        }
      }

    } ## Ends loop over j (columns).
  }   ## Ends loop over i (rows).

  if (use.legend) {
    plot(1, xlim = c(0, 1), ylim = c(0, 1), xlab = "", ylab = "",
         type = "n", axes = FALSE, xpd = NA)
    if (legend.location == "top") {
      xjust <- 0.5
      yjust <- 0
      legend.x <- .5
      legend.y <- 1
    } else {
      xjust <- 1
      yjust <- .5
      legend.x <- 1
      legend.y <- .5
    }

    if (style == "density") {
      legend(legend.x,
             legend.y,
             legend = names(list.of.arrays),
             col = color,
             lty = 1:length(list.of.arrays),
             bg = "white",
             horiz = (legend.location == "top"),
             xjust = xjust,
             yjust = yjust,
             cex = legend.cex,
             xpd = NA,
             bty = "n")
    } else if (style == "box") {
      legend(legend.x,
             legend.y,
             legend = names(list.of.arrays),
             fill = color,
             bg = "white",
             horiz = (legend.location == "top"),
             xjust = xjust,
             yjust = yjust,
             cex = legend.cex,
             xpd = NA,
             bty = "n")
    }
  }
  ## Add the main title, if one was supplied.
  if (main != "") {
    title(main = main, outer = TRUE, cex.main = 2, line = 5)
  }
}
