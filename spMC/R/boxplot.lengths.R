boxplot.lengths <-
function (x, ..., log = FALSE, zeros.rm = TRUE) {
  newlen <- x$length + 0
  if (zeros.rm & x$zeros) {
    newlen <- x$length + x$maxcens
  }
  if (log) newlen <- log(newlen)
  boxplot(newlen ~ x$categories, ...)
}

