#
#  Copyright (C) 2004-2012 Friedrich Leisch and Bettina Gruen
#  $Id: plot-flexmix.R 4922 2013-09-03 13:32:45Z gruen $
#

determine_y <- function(h, root) {
  y <- h$counts
  if (root) y <- sqrt(y)
  return(y)
}

panel.rootogram <- 
function (x, breaks, equal.widths = TRUE, nint = max(round(log2(length(x)) + 1), 3), alpha = plot.polygon$alpha,
          col = plot.polygon$col, border = plot.polygon$border, 
          lty = plot.polygon$lty, lwd = plot.polygon$lwd, subscripts, groups, mark, root = TRUE, markcol, ...) 
{
    x <- as.numeric(x)
    plot.polygon <- lattice::trellis.par.get("plot.polygon")
    grid::grid.lines(x = c(0.05, 0.95), y = grid::unit(c(0, 0), "native"), 
                     gp = grid::gpar(col = border, lty = lty, lwd = lwd, alpha = alpha),
               default.units = "npc")
    if (length(x) > 0) {
        if (is.null(breaks)) {
            breaks <- if (equal.widths) 
              lattice::do.breaks(range(x, finite = TRUE), nint)
            else quantile(x, 0:nint/nint, na.rm = TRUE)
        }
        h <- hist.constructor(x, breaks = breaks, plot = FALSE, ...)
        y <- determine_y(h, root)
        if (!is.null(mark)) {
          h1 <- hist.constructor(x[groups[subscripts] == mark], breaks = h$breaks, plot = FALSE, ...)
          y1 <- determine_y(h1, root)
        }
        nb <- length(breaks)
        if (length(y) != nb - 1) 
            warning("problem with hist computations")
        if (nb > 1) {
            lattice::panel.rect(x = breaks[-nb], y = 0, height = y, width = diff(breaks), 
                col = col, alpha = alpha, border = border, lty = lty, 
                lwd = lwd, just = c("left", "bottom"))
            if (!is.null(mark)) lattice::panel.rect(x = breaks[-nb], y = 0, height = y1, width = diff(breaks),
                                                    col = markcol, alpha = alpha, border = border, lty = lty, 
                                                    lwd = lwd, just = c("left", "bottom"))
        }
    }
}

prepanel.rootogram <- 
function (x, breaks, equal.widths = TRUE, nint = max(round(log2(length(x)) + 1), 3), root = TRUE, ...) 
{
  if (length(x) < 1) 
    list(xlim = NA, ylim = NA, dx = NA, dy = NA)
  else {
    if (is.factor(x)) {
      isFactor <- TRUE
      xlimits <- levels(x)
    }
    else isFactor <- FALSE
    if (!is.numeric(x)) 
      x <- as.numeric(x)
    if (is.null(breaks)) {
      breaks <- if (equal.widths) 
        lattice::do.breaks(range(x, finite = TRUE), nint)
      else quantile(x, 0:nint/nint, na.rm = TRUE)
    }
    h <- hist.constructor(x, breaks = breaks, plot = FALSE, ...)
    y <- determine_y(h, root)
    list(xlim = if (isFactor) xlimits else range(x, breaks, 
           finite = TRUE), ylim = range(0, y, finite = TRUE), 
         dx = 1, dy = 1)
  }
}


setMethod("plot", signature(x="flexmix", y="missing"),
function(x, y, mark=NULL, markcol=NULL, col=NULL, 
         eps=1e-4, root=TRUE, ylim=TRUE, main=NULL, xlab = "", ylab = "",
         as.table = TRUE, endpoints = c(-0.04, 1.04), ...){

    k <- length(x@prior)

    if(is.null(markcol)) markcol <- FullColors[5]
    if(is.null(col)) col <- LightColors[4]

    if(is.null(main)){
        main <- ifelse(root,
                      "Rootogram of posterior probabilities",
                      "Histogram of posterior probabilities")
        main <- paste(main, ">", eps)
    }
    groupfirst <- if (length(x@group)) !duplicated(x@group) else TRUE
    if (is.null(x@weights))
      z <- data.frame(posterior = as.vector(x@posterior$scaled[groupfirst,,drop=FALSE]),
                      component = factor(rep(seq_len(x@k), each = nrow(x@posterior$scaled[groupfirst,,drop=FALSE])),
                        levels = seq_len(x@k), labels = paste("Comp.", seq_len(x@k))),
                      cluster = rep(as.vector(x@cluster[groupfirst]), k))
    else 
      z <- data.frame(posterior = rep(as.vector(x@posterior$scaled[groupfirst,,drop=FALSE]),
                        rep(x@weights[groupfirst], k)),
                      component = factor(rep(seq_len(x@k), each = sum(x@weights[groupfirst])),
                        seq_len(x@k), paste("Comp.", seq_len(x@k))),
                      cluster = rep(rep(as.vector(x@cluster[groupfirst]), x@weights[groupfirst]), k))

    panel <- function(x, subscripts, groups, ...)
      panel.rootogram(x, root = root, mark = mark, col = col, markcol = markcol,
                      subscripts = subscripts, groups = groups, ...)
    prepanel <- function(x, ...) prepanel.rootogram(x, root = root, ...)
    z <- subset(z, posterior > eps)

    if (is.logical(ylim)) {
      scales <- if (ylim) list() else list(y = list(relation = "free"))
      hh <- lattice::histogram(~ posterior | component, data = z, main = main,  ylab = ylab, xlab = xlab, groups = cluster, 
                               panel = panel, prepanel = prepanel, scales = scales, as.table = as.table, endpoints = endpoints, ...)
    }
    else hh <- lattice::histogram(~ posterior | component, data = z, main = main, ylab = ylab, xlab = xlab, groups = cluster, 
                         ylim = ylim, panel = panel, prepanel = prepanel, as.table = as.table, endpoints = endpoints, ...)
    if (root) {
      hh$yscale.components <- function (lim, packet.number = 0, packet.list = NULL, right = TRUE, ...) 
        {
          comps <- calculateAxisComponents(lim, packet.list = packet.list, 
                                           packet.number = packet.number, ...)
          comps$at <- sqrt(seq(min(comps$at)^2, max(comps$at)^2, length.out = length(comps$at)))
          comps$labels <- format(comps$at^2, trim = TRUE)
          list(num.limit = comps$num.limit, left = list(ticks = list(at = comps$at, 
                                                            tck = 1), labels = list(at = comps$at, labels = comps$labels, 
                                                                          cex = 1, check.overlap = comps$check.overlap)), right = right)
        }
    }
    hh
  })
