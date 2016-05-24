summary.lengths <-
function (object, ..., zeros.rm = TRUE) {
  newlen <- object$length + 0
  if (zeros.rm & object$zeros) {
    newlen <- object$length + object$maxcens
  }
  res <- tapply(newlen, object$categories, function(p) summary(p, ...))
  class(res) <- "summary.lengths"
  return(res)
}

