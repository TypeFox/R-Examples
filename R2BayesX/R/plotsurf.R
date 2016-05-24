plotsurf <- function(data, x = 2, y = 3, z = 4, mode = 1, ...)
{
  if(is.character(data)) {
    stopifnot(file.exists(data <- path.expand(data)))
    sep <- list(...)$sep
    sep <- if(is.null(sep)) ""
    data <- read.table(data, header = TRUE, sep = sep)
  }
  data <- data[, c(x, y, z)]
  plot3d(x = data, c.select = 3:ncol(data),
    image = 2 %in% mode, contour = 3 %in% mode, ...)
}

