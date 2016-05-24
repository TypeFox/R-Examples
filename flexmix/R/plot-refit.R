#
#  Copyright (C) 2004-2012 Friedrich Leisch and Bettina Gruen
#  $Id: plot-refit.R 4859 2012-12-18 08:42:33Z gruen $
#

prepanel.default.coef <- function (x, y, subscripts, groups=NULL, horizontal = TRUE, nlevels, origin = NULL, 
                                   ...) 
{
  if (any(!is.na(x) & !is.na(y))) {
    if (horizontal) {
      if (!is.factor(y)) {
        if (missing(nlevels)) 
          nlevels <- length(unique(y))
        y <- factor(y, levels = seq_len(nlevels))
      }
      if (!is.null(groups)) {
        if (!is.numeric(x)) stop("x must be numeric")
              x <- rep(x, each = 2) + rep(groups[subscripts], each = 2) *c(-1,1)
      }
      list(xlim = if (is.numeric(x)) range(x, origin, finite = TRUE) else levels(x), 
           ylim = levels(y), yat = sort(unique(as.numeric(y))), 
           dx = 1, dy = 1)
    }
    else {
      if (!is.factor(x)) {
        if (missing(nlevels)) 
          nlevels <- length(unique(x))
        x <- factor(x, levels = seq_len(nlevels))
      }
      if (!is.null(groups)) {
        if (!is.numeric(y)) stop("y must be numeric")
        y <- rep(as.numeric(y), each = 2) + rep(groups[subscripts], each = 2) *c(-1,1)
      }
      list(xlim = levels(x), xat = sort(unique(as.numeric(x))), 
           ylim = if (is.numeric(y)) range(y, origin, finite = TRUE) else levels(y), 
           dx = 1, dy = 1)
    }
  }
  else list(xlim = c(NA, NA), ylim = c(NA, NA), dx = 1, dy = 1)
}

panel.coef <- function(x, y, subscripts, groups, significant = NULL, horizontal = TRUE, 
                       lwd = 2, col, col.line = c("black", "grey"), ...)
{
  col.sig <- rep(col.line[1], length(x))
  if (!is.null(significant)) {
    if (missing(col)) col <-   c("grey", "white")
    col.fill <- rep(col[1], length(x))
    col.sig[!significant[subscripts]] <- col.line[2]
    col.fill[!significant[subscripts]] <- col[2]
  }
  else if (missing(col)) col.fill <- "grey" else col.fill <- col
  lattice::panel.barchart(x, y, border = col.sig, col = col.fill, horizontal = horizontal, ...)
  if (!missing(groups)) {
    if (horizontal) {
      z <- x + rep(c(-1,1), each = length(x)) * matrix(rep(groups[subscripts], 2), ncol = 2)
      for (i in seq_along(x)) {
        lattice::panel.xyplot(z[i,], rep(y[i], 2), type = "l", col = col.sig[i], lwd = lwd)
      }
    }
    else {
      z <- y + rep(c(-1,1), each = length(y)) * matrix(rep(groups[subscripts], 2), ncol = 2)
      for (i in seq_along(y)) {
        lattice::panel.xyplot(rep(x[i], 2), z[i,], type = "l", col = col.sig[i], lwd = lwd)
      }
    }
  }
}

getCoefs <- function(x, alpha = 0.05, components, ...) {
  names(x) <- sapply(names(x), function(z) strsplit(z, "Comp.")[[1]][2])
  x <- x[names(x) %in% components]
  Comp <- lapply(names(x), function(n) 
                 data.frame(Value = x[[n]][,1],
                            SD = x[[n]][,2] * qnorm(1-alpha/2),
                            Variable = rownames(x[[n]]),
                            Component = n,
                            Significance = x[[n]][,4] <= alpha))
  do.call("rbind", Comp)
}

setMethod("plot", signature(x="FLXRoptim", y="missing"),
function(x, y, model = 1, which = c("model", "concomitant"),
         bycluster=TRUE, alpha=0.05, components, labels=NULL,
         significance = FALSE, xlab = NULL, ylab = NULL,
         ci = TRUE, scales = list(), as.table = TRUE, horizontal = TRUE, ...)
{
    which <- match.arg(which)
    if (missing(components)) components <- seq_len(x@k)
    plot.data <- if (which == "model") getCoefs(x@components[[model]], alpha, components) else getCoefs(x@concomitant, alpha, components)
    if (!is.null(labels)) plot.data$Variable <- factor(plot.data$Variable, labels = labels)
    plot.data$Component <- with(plot.data, factor(Component, sort(unique(Component)), labels = paste("Comp.", sort(unique(Component)))))
    if (bycluster) {
      formula <- if (horizontal) Variable ~ Value | Component else Value ~ Variable | Component
      plot.data$Variable <- with(plot.data, factor(Variable, levels = rev(unique(Variable))))
    }
    else {
      formula <- if (horizontal) Component ~ Value | Variable else Value ~ Component | Variable
      plot.data$Component <- with(plot.data, factor(Component, levels = rev(levels(Component))))
    }
    groups <- if (ci) plot.data$SD else NULL
    significant <- if (significance) plot.data$Significance else NULL
    lattice::xyplot(formula, data = plot.data, xlab = xlab, ylab = ylab, origin = 0, horizontal = horizontal,
                    scales = scales, as.table = as.table, significant = significant,
                    groups = groups, prepanel = function(...) prepanel.default.coef(...),
                    panel = function(x, y, subscripts, groups, ...)
                    panel.coef(x, y, subscripts, groups, ...), ...)
})

