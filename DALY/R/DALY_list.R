## Create 'DALY_list' objects

DALY_list <-
function(...){
  ## catch arguments
  x <- list(...)

  ## evaluate arguments
  if (length(x) == 1)
    stop("At least two 'DALY' objects must be provided")
  if (!all(unlist(lapply(x, function(x) class(x)[1] == "DALY")))) 
    stop("All arguments must be 'DALY' objects")

  ## define S3 class
  class(x) <- "DALY_list"

  ## return 'DALY_list'
  return(x)
}

print.DALY_list <-
function(x, ...){
  for (i in seq_along(x))
    print(x[[i]], ...)
}

## PLOT method for class 'DALY_list'
## generates stacked aggregated YLL/YLD barplot
## with overall DALY credibility interval
plot.DALY_list <-
function(x, prob = 0.95, sort = TRUE, names = NULL,
         bars = TRUE, col = c("grey90", "white"),
         error_bars = TRUE, eb_col = "black",
         grid = TRUE, ...){

  ## check input values
  if (!bars & !error_bars)
    stop("'bars' and 'error_bars' cannot both be FALSE")

  ## number of scenarios
  n <- length(x)

  ## check names
  if (is.null(names)){
    names <- sapply(x, function(x) x$name)
  } else {
    if (length(names) != n)
      stop("Length of 'names' does not match number of DALY objects (",
           n, ")")
  }

  ## obtain summary statistics
  y <- lapply(x, aggregate)

  ## generate plot
  DALY_barplot(y, n, prob, sort, names,
               bars, col, error_bars, eb_col, grid, ...)
}