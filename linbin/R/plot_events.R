#' Plot Events as Bar Plots
#' 
#' Plots an event table as a grid of bar plots.
#' 
#' Given a groupping variable for the rows of the event table (e.g., groups of bins of different sizes used in \code{\link{sample_events}}), and groups of columns to plot, bar plots are drawn in a grid for each combination of event and column groups. In each plot, the specified event table columns are plotted together as stacked bars. Negative and positive values are stacked seperately from the \code{y = 0} baseline. Events with \code{NA} are not shown, differentiating them from zero-valued events which are drawn as thin black lines. Point events are drawn as thin vertical lines. Overlapping events are drawn as overlapping bars, so it is best to use \code{\link{sample_events}} with non-overlapping bins to flatten the data before plotting.
#'
#' @param e An event table.
#' @param group.col Name or index of column defining the event groupping for plotting. If \code{NULL}, the events are treated as one group. Group \code{NA} is not plotted.
#' @param groups Vector of values from \code{group.col} specifying which groups to plot. If \code{NULL}, all groups are plotted by order of first appearance in \code{group.col}.
#' @param data.cols Names or indices of columns to plot, given as a list of character or numeric vectors. If multiple columns are specified, their bars are stacked together in one plot. Names are interpreted as regular expressions (\code{\link{regex}}) matching full column names. If \code{NULL}, all columns not named \code{from}, \code{to}, or \code{group.col} are each plotted individually in order of appearance.
#' @param dim The row and column dimensions of the grid. If \code{NULL}, the grid is column groups (rows) by event groups (columns) if \code{byrow = TRUE}, and event groups (rows) by column groups (columns) if \code{byrow = FALSE}.
#' @param byrow Plots are added by column group, then bin group. If \code{TRUE}, plots are added by rows, rather than columns, to the grid.
#' @param main Titles for each plot. If \code{NULL}, plots are titled by the column names, pasted together with seperator " + ". Set \code{main = NA} to not title the plots.
#' @param xlabs,ylabs Labels arranged at equal intervals along the bottom and left side of the plot grid. These are drawn in the outer margins of the figure, so \code{oma[1]} and \code{oma[2]} must be non-zero.
#' @param xlim,ylim Limits for the x and y axes of all plots. If \code{NULL}, limits are set to the range of the data and the y limits extended as needed to include 0.
#' @param xticks,yticks The positions of x and y tick marks for all plots. If \code{NULL}, only the min and max x and y are ticked (and 0 as needed for y). If \code{\link{axTicks}}, that function will be used to calculate R default tick mark positions. If \code{NA}, no ticks are drawn.
#' @param xtick.labels,ytick.labels The labels for the x and y tick marks, coerced to character vectors and recyled as necessary. If \code{NULL}, the positions of the ticks are used as the labels, formatted with \code{sigfigs}. If \code{NA}, the tick marks are not labeled.
#' @param plot.grid If \code{TRUE}, a lined horizontal grid is plotted at the yticks.
#' @param sigfigs The maximum significant figures of the x and y axis labels.
#' @param col Color(s) for the bars in each plot. If \code{NA}, bars are transparent. If \code{NULL}, a grey palette is used.
#' @param border Color(s) for bar borders in each plot. If \code{NA}, borders are omitted.
#' @param lty Line type(s) for bar borders in each plot.
#' @param lwd Line width(s) for bar borders in each plot.
#' @param xpd Logical value or \code{NA}. If \code{FALSE}, all plotting is clipped to the plot region, if \code{TRUE}, all plotting is clipped to the figure region, and if \code{NA}, all plotting is clipped to the device region.
#' @param mar Numerical vector of the form c(bottom, left, top, right) giving the size of the inner margins of each plot in lines of text.
#' @param oma Numeric vector of the form c(bottom, left, top, right) giving the size of the outer figure margins in lines of text.
#' @param ... Additional arguments passed to \code{\link{plot}}.
#' @seealso \code{\link{seq_events}} for generating groups of sequential bins, \code{\link{sample_events}} to populate groups of bins with event data.
#' @export
#' @import graphics
#' @examples
#' e <- events(from = c(0, 10, 15, 25), to = c(10, 20, 25, 40), length = c(10, 10, 10, 15),
#'             x = c(1, 2, 1, 1), f = c('a', 'b', 'a', 'a'))
#' bins <- seq_events(event_coverage(e), c(8, 4, 2, 1))
#' e.bins <- sample_events(e, bins, list(sum, c('x', 'length')), scaled.cols = 'length')
#' plot_events(e.bins, group.col = 'group')
plot_events <- function(e, group.col = NULL, groups = NULL, data.cols = NULL, dim = NULL, byrow = TRUE, main = NULL, xlabs = character(), ylabs = character(), xlim = NULL, ylim = NULL, xticks = NULL, yticks = NULL, xtick.labels = NULL, ytick.labels = NULL, plot.grid = FALSE, sigfigs = c(3, 3), col = NULL, border = par("fg"), lty = par("lty"), lwd = par("lwd"), xpd = FALSE, mar = c(2.1, 2.75, 1.5, 0.5),  oma = c(2, 2, 2, 2), ...) {
  
  ### Initialize
  if (is.null(group.col)) {
    # all events in one group
    group.seq <- numeric(nrow(e))
    groups <- 0
  } else {
    group.seq <- e[[group.col]]
    if (is.null(groups)) {
      # all event groups
      groups <- unique(group.seq)
      groups <- groups[!is.na(groups)]
    } else {
      # user-specified event groups
      groups <- unique(groups)
    }
  }
  if (is.null(data.cols)) {
    group.col <- if_else(is.character(group.col), group.col, names(e)[group.col])
    data.cols <- c(list(), setdiff(names(e), c("from", "to", group.col)))
  } else {
    # FIXME: Make new helper function for retrieving column names or indices from user input?
    # Regex to indices
    data.cols <- rapply(as.list(data.cols), rgrep_exact, x = names(e), classes = "character", how = "replace")
    # FIXME: Why is the unlisting needed?
    data.cols <- lapply(data.cols, unlist)
    # Indices to names
    ind2name <- function(x) names(e)[x]
    data.cols <- rapply(data.cols, ind2name, classes = "ANY", how = "replace")
    data.cols <- lapply(data.cols, unlist)
  }
  
  # grid dimensions
  nr <- length(data.cols)
  nc <- length(groups)
  n <- nr * nc
  if (byrow) {
    par(mfrow = if_else(is.null(dim), c(nr, nc), dim))
  } else {
    par(mfcol = if_else(is.null(dim), c(nc, nr), dim))
  }
  par(mar = mar, oma = oma)
  main <- if_else(is.null(main), rep(sapply(data.cols, paste, collapse = " + "), rep(nc, nr)), main)
  
  ### Plot each combination
  j <- 1
  for (i in seq_along(data.cols)) {
    plot.col = if_else(is.null(col), grDevices::grey.colors(length(data.cols[[i]])), col)
    for (group in groups) {
      ind <- group == group.seq
      plot_events_single(e[ind, , drop = FALSE], data.cols[[i]], main = main[j], xlim = xlim, ylim = ylim, xticks = xticks, yticks = yticks, sigfigs = sigfigs, col = plot.col, border = border, lty = lty, lwd = lwd, xpd = xpd, plot.grid = plot.grid, xtick.labels = xtick.labels, ytick.labels = ytick.labels, ...)
      j <- j + 1
    }
  }
  
  ### Outer axis labels
  if (!is.na(xlabs) && length(xlabs)) {
    n = length(xlabs)
    left = (0:(n - 1) / n) + par("plt")[1] * diff(par("fig")[1:2])
    right = ((1:n) / n) - (1 - par("plt")[2]) * diff(par("fig")[1:2])
    mtext(xlabs, side = 1, outer = TRUE, at = rowMeans(cbind(left, right)), cex = 0.95, line = max(0, par("oma")[1] - 2))
  }
  if (!is.na(ylabs) && length(ylabs)) {
    n = length(ylabs)
    bottom = (0:(n - 1) / n) + par("plt")[3] * diff(par("fig")[3:4])
    top = ((1:n) / n) - (1 - par("plt")[4]) * diff(par("fig")[3:4])
    mtext(ylabs, side = 2, outer = TRUE, at = rowMeans(cbind(bottom, top)), cex = 0.95, line = max(0, par("oma")[2] - 2))  
  }
}

#' Plot Events as Bars
#' 
#' Plots event table columns as vertical bars.
#' 
#' The specified event table columns are plotted together as stacked bars. Negative and positive values are stacked seperately from the \code{y = 0} baseline. Events with \code{NA} are not shown, differentiating them from zero-valued events which are drawn as thin black lines. Point events are drawn as thin vertical lines. Overlapping events are drawn as overlapping bars, so it is best to use \code{\link{sample_events}} with non-overlapping bins to flatten the data before plotting.
#' 
#' @param e An event table.
#' @param cols Names or indices of the event table columns to plot together as stacked bars.
#' @param xlim,ylim Limits for the x and y axes. If \code{NULL}, limits are set to the range of the data and the y limits extended as needed to include 0.
#' @param xticks,yticks The values to label on the x and y axes. If \code{NULL}, only the min and max x and y are labeled (and 0 as needed for y). If \code{\link{axTicks}}, the function will be used to generate R default tick marks.
#' @param xtick.labels,ytick.labels Labels for the x and y tick positions.
#' @param main An overall title for the plot.
#' @param xlab,ylab Titles for the x and y axes.
#' @param plot.grid If \code{TRUE}, a lined horizontal grid is plotted at the yticks.
#' @param sigfigs The maximum significant figures to use for the x and y axis labels.
#' @param col Color(s) for the bars. If \code{NULL}, bars are transparent. By default, a grey palette is used.
#' @param border Color(s) for bar borders. Use border = NA to omit borders.
#' @param lty Line type(s) for bar borders.
#' @param lwd Line width(s) for bar borders.
#' @param xpd Logical value or \code{NA}. If \code{FALSE}, all plotting is clipped to the plot region, if \code{TRUE}, all plotting is clipped to the figure region, and if \code{NA}, all plotting is clipped to the device region.
#' @param ... Additional arguments passed to \code{\link{plot}}.
#' @seealso \code{\link{plot_events}}.
#' @keywords internal
#' @import graphics
plot_events_single <- function(e, cols, xlim = NULL, ylim = NULL, xticks = NULL, yticks = NULL, xtick.labels = NULL, ytick.labels = NULL, main = NA, xlab = NA, ylab = NA, plot.grid = FALSE, sigfigs = c(3, 3), col = grDevices::grey.colors(length(cols)), border = par("fg"), lty = par("lty"), lwd = par("lwd"), xpd = FALSE, ...) {
  
  # Compute plot limits
  if (is.null(xlim)) {
    xlim <- c(min(e$from), max(e$to))
  }
  if (is.null(ylim)) {
    sum.neg <- function(x, ...) sum(x[x < 0], ...)
    sum.pos <- function(x, ...) sum(x[x > 0], ...)
    y.neg <- apply(e[cols], 1, sum.neg, na.rm = TRUE)
    y.pos <- apply(e[cols], 1, sum.pos, na.rm = TRUE)
    ylim <- range(c(0, y.neg, y.pos))
  }
  
  # Initialize plot
  plot(xlim, ylim, type = 'n', axes = FALSE, xlim = xlim, ylim = ylim, main = main, ylab = ylab, xlab = xlab, ...)
  
  # Axis ticks
  if (is.null(xticks)) {
    xticks <- xlim
  } else if (identical(xticks, axTicks)) {
    xticks <- axTicks(1)
  }
  if (is.null(yticks)) {
    yticks <- unique(c(0, ylim))
  } else if (identical(yticks, axTicks)) {
    yticks <- unique(c(0, axTicks(2)))
  }
  if (length(sigfigs) == 1) {
    sigfigs <- rep(sigfigs, 2)
  }
  if (is.null(xtick.labels)) {
    if (is.null(sigfigs)) {
      xtick.labels <- xticks
    } else {
      xtick.labels <- prettyNum(signif(xticks, sigfigs[1]), big.mark = ",")
    }
  } else {
    xtick.labels <- rep_len(xtick.labels, length(xticks))
  }
  if (is.null(ytick.labels)) {
    if (is.null(sigfigs)) {
      ytick.labels <- yticks
    } else {
      ytick.labels <- prettyNum(signif(yticks, sigfigs[2]), big.mark = ",")
    }
  } else {
    ytick.labels <- rep_len(ytick.labels, length(yticks))
  }
  
  # Draw axes, ticks, and labels
  if (!identical(NA, xticks)) {
    axis(side = 1, at = xticks, labels = xtick.labels, col = NA, col.ticks = par("fg"))
  }
  if (!identical(NA, yticks)) {
    if (plot.grid) {
      axis(side = 2, at = yticks, labels = FALSE, tck = 1, lty = 3)
    }
    axis(side = 2, at = yticks, labels = ytick.labels, las = 1)
  }
  
  # Plot stacked barplots
  # (for each set of bars, stack positive on positive and negative on negative)
  nc <- length(cols)
  ne <- nrow(e)
  col <- if_else(is.null(col), col, rep_len(col, nc))
  border <- if_else(is.null(border), border, rep_len(border, nc))
  lwd <- if_else(is.null(lwd), lwd, rep_len(lwd, nc))
  lty <- if_else(is.null(lty), lty, rep_len(lty, nc))
  base.pos <- rep(0, ne)
  base.neg <- rep(0, ne)
  for (i in seq_len(nc)) {
    h <- e[[cols[i]]]
    y0 <- ifelse(h > 0, base.pos, base.neg)
    y1 <- y0 + h
    rect(e$from, y0, e$to, y1, col = col[i], border = border[i], lwd = lwd[i], lty = lty[i], xpd = xpd)
    base.pos <- ifelse(!is.na(h) & h > 0, base.pos + h, 0)
    base.neg <- ifelse(!is.na(h) & h < 0, base.neg + h, 0)
  }
}