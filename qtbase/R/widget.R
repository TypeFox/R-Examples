### Conveniences for widgets

print.QWidget <- function(x, ...)
{
  x$show()
  NextMethod()
  invisible(x)
}
