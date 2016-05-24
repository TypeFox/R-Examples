iter.connection <- function(obj, ...) {
  s <- summary(obj)

  if (s$opened != "opened")
    stop("connection not opened")

  if (s$`can read` != "yes")
    stop("connection not readable")

  if (s$text == "binary")
    ireadBin(obj, ...)
  else
    ireadLines(obj, ...)
}
