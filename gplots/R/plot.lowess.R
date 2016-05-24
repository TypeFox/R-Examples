plotLowess <- function (formula, data=parent.frame(), ..., subset=parent.frame(),
                        col.lowess="red",
                        lty.lowess=2  )
{
  m <- match.call(expand.dots=TRUE)
  m[[1]] <- as.name("plot")
  eframe <- parent.frame()
  eval(m, eframe)
  m[[1]] <- as.name("lowess")
  lw <- eval(m, eframe)
  lines(lw, col=col.lowess, lty=lty.lowess)
  grid()

  invisible(lw)
}

plot.lowess <- function(x, y, ..., col.lowess="red", lty.lowess=2)
{
  m <- x$call
  m[[1]] <- quote(plot)
  m <- as.call(append(as.list(m), list(...)))
  eval(m,  envir = parent.frame())

  lines(x$x, x$y, col=col.lowess, lty=lty.lowess)
  grid()

  invisible(x)
}
