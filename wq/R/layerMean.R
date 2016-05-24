layerMean <-
function(d) {
    
  # Trapezoidal mean of scalar x versus z
  trapMean <- function(z, x) {
    # Handle NAs
    w <- na.omit(cbind(z, x))
    n <- nrow(w)
    if (identical(n, 0L))
      return(NA)
    z <- w[, 1]
    x <- w[, -1]
    z1 <- diff(z)
    x1 <- 0.5 * (x[-1] + x[-n])
    sum(z1 * x1)/(z[n] - z[1])
  }
  
  # Trapezoidal mean of vector d[, -1] vs d[, 1]
  # Handle single observations
  n <- nrow(d)
  if (is.null(n))
    return(d[-1])
  if (identical(n, 1L))
    return(as.numeric(d[, -1]))
  # Handle duplicates
  d <- aggregate(d[, -1], by = list(z = d[, 1]), mean, na.rm = TRUE)

  apply(d[, -1, drop = FALSE], 2, function(x) trapMean(z = d[, 1], x))
}
