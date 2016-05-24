`summary.growth` <-
function(object, ...) {
  if (!is(object, "growth")) stop("argument should be a growth object")
  print(object)
  invisible(object)
}

