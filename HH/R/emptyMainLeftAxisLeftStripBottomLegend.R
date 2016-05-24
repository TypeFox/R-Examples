emptyMainLeftAxisLeftStripBottomLegend <- function(x) {
  ## main title
  x <- update(x, main=" ")
  ## subtitle
  if (!is.null(x$sub)) x <- update(x, sub=" ")
  ## left tick labels
  if (is.list(x$y.limits))
    x$y.limits <- lapply(x$y.limits,
                         function(x) {
                           x[] <- ""
                           x
                         }
                         )
  else
    x$y.limits[] <- " "
  ## left strip
  x$par.strip.text$lines <- 0
  ## bottom legend
  x$legend$bottom$args$text[] <- " "
  x$legend$bottom$args$rect[] <- 0
  if (!is.null(x$legend$bottom$args$title))
    x$legend$bottom$args$title <- " "
  x
}

## CountRev <- emptyMainLeftAxisLeftStripBottomLegend(CountsPlot)

## source("c:/HOME/rmh/HH-R.package/HH/R/emptyMainLeftAxisLeftStripBottomLegend.R")
