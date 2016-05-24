## PLOT method for class 'DALY'
## generates stacked YLL/YLD barplot
## with overall DALY credibility interval

plot.DALY <-
function(x, prob = 0.95, sort = TRUE, names = NULL,
         bars = TRUE, col = c("grey90", "white"),
         error_bars = TRUE, eb_col = "black",
         grid = TRUE, ...){

  ## check input values
  if (!bars & !error_bars)
    stop("'bars' and 'error_bars' cannot both be FALSE")

  ## number of outcomes
  n <- length(x) - 2

  ## check names
  if (is.null(names)){
    names <- sapply(x, function(x) if(is.list(x)) x$name)[seq(n)]
  } else {
    if (length(names) != n)
      stop("Length of 'names' does not match number of outcomes (", n, ")")
  }

  ## obtain summary statistics
  y <- aggregate(x, by = "outcome")

  ## generate plot
  DALY_barplot(y, n, prob, sort, names,
               bars, col, error_bars, eb_col, grid, ...)
}