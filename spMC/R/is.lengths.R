is.lengths <-
function (object) {
  if (!is(object, "lengths")) return(FALSE)
  if (length(object) != 5) return(FALSE)
  if (!prod(c("categories", "direction", "maxcens", "length", "zeros") %in% names(object))) return(FALSE)
  if (!is.factor(object$categories)) return(FALSE)
  if (!is.numeric(object$direction)) return(FALSE)
  if (!is.numeric(object$maxcens)) return(FALSE)
  if (!is.numeric(object$length)) return(FALSE)
  if (!is.logical(object$zeros)) return(FALSE)
  return(TRUE)
}

