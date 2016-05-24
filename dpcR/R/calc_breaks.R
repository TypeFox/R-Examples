# Function for all cases when we need breaks calculated
calc_breaks <- function(vals, breaks = "Sturges", threshold = NULL) {
  if (!is.vector(vals))
    vals <- as.vector(sapply(vals, function(x) as.vector(x)))
  if (!is.null(threshold)) {
    min_v <- min(vals)
    vals <- vals[vals > threshold]
  }
  br <- hist(vals, breaks, plot = FALSE)[["breaks"]]
  if (!is.null(threshold))
    br <- c(min_v, br)
  br
}