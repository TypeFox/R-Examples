# mtplot.r
# Time-stamp: c:/x/rpack/kw/R/fcdf.r

mountainplot <- function (x, data, ...)
  UseMethod("mountainplot")

mountainplotyscale.components <- function(...) {
  ans <- yscale.components.default(...)
  ans$right <- ans$left
  foo <- ans$right$labels$at
  ans$right$labels$labels <- as.character(1-foo)
  ans
}

mountainplot.formula <- function (x, data = NULL,
                                  prepanel = "prepanel.mountainplot",
                                  panel = "panel.mountainplot",
                                  ylab = gettext("Folded Empirical CDF"),
                                  yscale.components = mountainplotyscale.components,
                                  scales = list(y = list(alternating = 3)),
                                  ...) {
  ccall <- match.call()
  ocall <- sys.call(sys.parent())
  ocall[[1]] <- quote(mountainplot)
  ccall$data <- data
  ccall$prepanel <- prepanel
  ccall$panel <- panel
  ccall$ylab <- ylab
  ccall$yscale.components <- yscale.components
  ccall$scales <- scales
  ccall[[1]] <- quote(lattice::densityplot)  # Why...?
  ans <- eval.parent(ccall)
  ans$call <- ocall
  ans
}

mountainplot.numeric <- function (x, data = NULL, xlab = deparse(substitute(x)), ...) {
    ccall <- match.call()
    ocall <- sys.call(sys.parent())
    ocall[[1]] <- quote(mountainplot)
    if (!is.null(ccall$data))
        warning("explicit 'data' specification ignored")
    ccall$data <- list(x = x)
    ccall$xlab <- xlab
    ccall$x <- ~x
    ccall[[1]] <- quote(mountainplot)  # See note from Felix Andrews
    ans <- eval.parent(ccall)
    ans$call <- ocall
    ans
}

prepanel.mountainplot <- function (x, ...) {
  ans <- prepanel.default.qqmath(x, distribution = qunif)
  with(ans, list(xlim = ylim, ylim = c(0, .5), dx = dy, dy = dx))
}

panel.mountainplot <- function (x, type = "s",
                                groups = NULL,
                                ref = TRUE, ...) {
  reference.line <- trellis.par.get("reference.line")
  if (ref) {
    reference.line <- trellis.par.get("reference.line")
    do.call(panel.abline, c(list(h = c(0, 1)), reference.line))
  }
  x <- as.numeric(x)
  distribution <- qunif
  nobs <- sum(!is.na(x))
  if (!is.null(groups)) {
    panel.superpose(x, y = NULL, type = type,
                    distribution = distribution, groups= groups,
                    panel.groups = panel.mountainplot, ...)
  }
  else if (nobs) {
    ypos <- seq_len(nobs)/(nobs+1)
    ypos <- ifelse(ypos<=.5, ypos, 1-ypos)
    panel.xyplot(x = sort(x), y = ypos, type = type, ...)
    panel.abline(h = c(.1,.25), col=reference.line$col, lty=2)
  }
}
