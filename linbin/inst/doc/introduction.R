## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE, 
  comment = '>',
  fig.align = 'center',
  fig.show = 'hold'
)

## ---- echo=FALSE---------------------------------------------------------
library(linbin)

## ------------------------------------------------------------------------
e <- simple

## ---- echo = FALSE, results = 'asis'-------------------------------------
knitr::kable(e)

## ---- results = 'asis'---------------------------------------------------
bins <- seq_events(event_range(e), length.out = 5)
e.bins <- sample_events(e, bins, list(mean, "x"), list(mean, "y", by = "factor", na.rm = TRUE))

## ---- echo = FALSE, results = 'asis'-------------------------------------
row.names(e.bins) = NULL
knitr::kable(e.bins)

## ---- fig.width = 4, fig.height = 4--------------------------------------
plot_events(e.bins, xticks = axTicks, border = par("bg"))

## ------------------------------------------------------------------------
events(from = c(0, 15, 25), to = c(10, 30, 35), x = 1, y = c('a', 'b', 'c'))

## ------------------------------------------------------------------------
as_events(1:3) # vector
as_events(cbind(1:3, 2:4)) # matrix
as_events(data.frame(start = 1:3, x = 1, stop = 2:4), "start", "stop") # data.frame

## ---- echo = FALSE, fig.width = 4, fig.height = 4------------------------
metrics <- rbind(
  cbind(event_range(e), n = 1, g = 1),
  cbind(event_coverage(e, closed = FALSE), n = 1, g = 2),
  cbind(event_overlaps(e), g = 4),
  cbind(event_gaps(e), n = 1, g = 3))
plot_events(metrics, group.col = "g", data.cols = "n", dim = c(4, 1), xlim = c(0, 100), 
            col = 'grey', main = c("Range", "Coverage", "Overlaps", "Gaps"), 
            lwd = 0.5, mar = rep(1.5, 4), oma = rep(1, 4))

## ------------------------------------------------------------------------
e <- elwha

## ------------------------------------------------------------------------
bins <- event_overlaps(e)

## ---- fig.width = 4.5, fig.height = 1.25---------------------------------
e.bins <- sample_events(e, bins, list(weighted.mean, "mean.width", "unit.length"), 
                       scaled.cols = "unit.length")
plot_events(e.bins, data.cols = "mean.width", col = "grey", border = "#666666", 
            ylim = c(0, 56), main = "", oma = rep(0, 4), mar = rep(0, 4), 
            xticks = NA, yticks = NA)

## ------------------------------------------------------------------------
bins <- seq_events(event_range(e), length.out = 33)

## ---- echo = FALSE, fig.width = 4.5, fig.height = 1.25-------------------
e.bins <- sample_events(e, bins, list(weighted.mean, "mean.width", "unit.length"), scaled.cols = "unit.length")
plot_events(e.bins, data.cols = "mean.width", col = "grey", border = "#666666", 
            ylim = c(0, 56), main = "", oma = rep(0, 4), mar = rep(0, 4), 
            xticks = NA, yticks = NA)

## ------------------------------------------------------------------------
bins <- seq_events(event_coverage(e), length.out = 20)

## ---- echo = FALSE, fig.width = 4.5, fig.height = 1.25-------------------
e.bins <- sample_events(e, bins, list(weighted.mean, "mean.width", "unit.length"), scaled.cols = "unit.length")
plot_events(e.bins, data.cols = "mean.width", col = "grey", border = "#666666", 
            ylim = c(0, 56), main = "", oma = rep(0, 4), mar = rep(0, 4), 
            xticks = NA, yticks = NA)

## ------------------------------------------------------------------------
e.filled <- fill_event_gaps(e, max.length = 1) # fill small gaps first
bins <- seq_events(event_coverage(e.filled), length.out = 20, adaptive = TRUE)

## ---- echo = FALSE, fig.width = 4.5, fig.height = 1.25-------------------
e.bins <- sample_events(e, bins, list(weighted.mean, "mean.width", "unit.length"), scaled.cols = "unit.length")
plot_events(e.bins, data.cols = "mean.width", col = "grey", border = "#666666", 
            ylim = c(0, 56), main = "", oma = rep(0, 4), mar = rep(0, 4), 
            xticks = NA, yticks = NA)

## ------------------------------------------------------------------------
e <- simple
bins <- seq_events(event_range(e), length.out = 1)

## ------------------------------------------------------------------------
e.bins <- sample_events(e, bins, list(sum, c('x', 'y'), na.rm = TRUE), scaled.cols = c('x', 'y'))

## ---- echo = FALSE, results = 'asis'-------------------------------------
row.names(e.bins) = NULL
knitr::kable(e.bins)

## ------------------------------------------------------------------------
e.bins <- sample_events(e, bins, list(weighted.mean, 'x', 'y', na.rm = TRUE))

## ---- echo = FALSE, results = 'asis'-------------------------------------
row.names(e.bins) = NULL
knitr::kable(e.bins)

## ------------------------------------------------------------------------
fun <- function(x) paste0(unique(x), collapse = '.')
e.bins <- sample_events(e, bins, list(fun, 'factor'))

## ---- echo = FALSE, results = 'asis'-------------------------------------
row.names(e.bins) = NULL
knitr::kable(e.bins)

## ----fig.height = 4, fig.width = 5---------------------------------------
e <- simple
bins <- seq_events(event_range(e), length.out = c(16, 4, 2)) # appends a "group" column
e.bins <- sample_events(e, bins, list(sum, c('x', 'y'), na.rm = TRUE))
plot_events(e.bins, group.col = 'group')

