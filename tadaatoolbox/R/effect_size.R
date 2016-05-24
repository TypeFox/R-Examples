#' Simple Effect Size Calculation for t-Tests
#'
#' @param data A \code{data.frame}.
#' @param response The response variable (dependent).
#' @param group The group variable, usually a \code{factor}.
#' @param na.rm If \code{TRUE} (default), missing values are dropped.
#' @return \code{numeric} of length 1.
#' @export
#' @import stats
#' @examples
#' df <- data.frame(x = runif(100), y = sample(c("A", "B"), 100, TRUE))
#' effect_size_t(df, "x", "y")
effect_size_t <- function(data, response, group, na.rm = TRUE){
  # Check the type of the group
  if (is.factor(data[[group]])) {
    groups <- levels(data[[group]])
  } else {
    groups <- unique(data[[group]])
  }

  # Subset groups of response
  x <- data[data[[group]] == groups[1], ][[response]]
  y <- data[data[[group]] == groups[2], ][[response]]

  # Kick out NAs if specified
  if (na.rm) {
    x <- x[!is.na(x)]
    y <- y[!is.na(y)]
  }

  # Get stats for each group
  n1   <- length(x)
  var1 <- var(x)

  n2   <- length(y)
  var2 <- var(y)

  # Calculate pooled variance and difference of means
  s    <- sqrt(sum(n1 * var1, n2 * var2)/((n1 + n2) - 2))
  m_d  <- mean(x) - mean(y)

  # Get effect size (absolute)
  d    <- abs(m_d/s)

  return(d)
}
