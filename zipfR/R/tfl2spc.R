tfl2spc <- function (tfl)
{
  if (! inherits(tfl, "tfl")) stop("argument must be object of class 'tfl'")
  if (attr(tfl, "incomplete")) stop("incomplete type frequency lists are not supported")

  f <- tfl$f
  if (!is.integer(f)) {
    if (! all(tfl$f == floor(tfl$f))) stop("type frequencies in 'tfl' must be integer values")
    f <- as.integer(f)
  }

  x <- rle(sort(tfl$f))
  spc(Vm=x$lengths, m=x$values)
}
