plotnonp <- function(data, x = 2, y = c(3, 4, 5, 7, 8),
  lty = c(1, 3, 2, 2, 3), cols = NULL, month = NULL, year = NULL,
  step = 12, ...)
{
  plot2d(x = data, c.select = c(x, y), col.lines = cols, lty = lty,
    month = month, year = year, step = step, ...)
}
