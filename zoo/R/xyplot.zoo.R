panel.plot.default <- function(x, y, subscripts, groups, panel = panel.xyplot,
  col = 1, type = "p", pch = 20, lty = 1, lwd = 1, ...)
{
  col <- rep(as.list(col), length.out = nlevels(groups))
  type <- rep(as.list(type), length.out = nlevels(groups))
  pch <- rep(as.list(pch), length.out = nlevels(groups))
  lty <- rep(as.list(lty), length.out = nlevels(groups))
  lwd <- rep(as.list(lwd), length.out = nlevels(groups))

  for(g in 1:nlevels(groups)) {
    idx <- g == unclass(groups[subscripts])
    if (any(idx)) panel(x[idx], y[idx], ...,
      col = col[[g]], type = type[[g]], pch = pch[[g]],
      lty = lty[[g]], lwd = lwd[[g]])
  }
  .Deprecated(msg="panel.plot.default is no longer needed, just use panel.xyplot etc")
}

panel.plot.custom <- function(...) {
  args <- list(...)
  function(...) {
    dots <- list(...)
    do.call("panel.plot.default", modifyList(dots, args))
  }
}

xyplot.its <-
xyplot.ts <-
xyplot.zoo <- function(x, data, ...)
{
    obj <- lattice::xyplot.ts(as.zoo(x), ...)
    obj$call <- match.call()
    obj
}

xyplot.tis <- function(x, data, ...)
{
    x <- aggregate(as.zoo(x), tis::POSIXct, identity)
    obj <- lattice::xyplot.ts(x, ...)
    obj$call <- match.call()
    obj
}


llines.its <-
llines.tis <-
llines.zoo <- function(x, y = NULL, ...)
{
    if (!is.null(y)) {
        llines(coredata(x), y = y, ...)
    } else {
        llines(coredata(time(x)), y = coredata(x), ...)
    }
}

lpoints.its <-
lpoints.tis <-
lpoints.zoo <- function(x, y = NULL, ...)
{
    if (!is.null(y)) {
        lpoints(coredata(x), y = y, ...)
    } else {
        lpoints(coredata(time(x)), y = coredata(x), ...)
    }
}

ltext.its <-
ltext.tis <-
ltext.zoo <- function(x, y = NULL, ...)
{
    if (!is.null(y)) {
        ltext(coredata(x), y = y, ...)
    } else {
        ltext(coredata(time(x)), y = coredata(x), ...)
    }
}


panel.lines.ts <- 
panel.lines.its <-
panel.lines.tis <-
panel.lines.zoo <- function(x, ...) {
  x <- as.zoo(x)
  panel.lines(time(x), coredata(x), ...)
  .Deprecated("panel.lines")
}

panel.points.ts <- 
panel.points.its <-
panel.points.tis <-
panel.points.zoo <- function(x, ...) {
  x <- as.zoo(x)
  panel.points(time(x), coredata(x), ...)
  .Deprecated("panel.points")
}

panel.text.ts <- 
panel.text.its <-
panel.text.tis <-
panel.text.zoo <- function(x, ...) {
  x <- as.zoo(x)
  panel.text(time(x), coredata(x), ...)
  .Deprecated("panel.text")
}

panel.segments.ts <- 
panel.segments.its <-
panel.segments.tis <-
panel.segments.zoo <- function(x0, x1, ...) {
  x0 <- as.zoo(x0)
  x1 <- as.zoo(x1)
  panel.segments(time(x0), coredata(x0), time(x1), coredata(x1), ...)
}

panel.rect.ts <- 
panel.rect.its <-
panel.rect.tis <-
panel.rect.zoo <- function(x0, x1, ...) {
  x0 <- as.zoo(x0)
  x1 <- as.zoo(x1)
  panel.rect(time(x0), coredata(x0), time(x1), coredata(x1), ...)
}

panel.polygon.ts <- 
panel.polygon.its <-
panel.polygon.tis <-
panel.polygon.zoo <- function(x, ...) {
  x <- as.zoo(x)
  panel.polygon(time(x), coredata(x), ...)
}
