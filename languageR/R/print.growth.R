`print.growth` <-
function(x, ...) {
  if (!is(x, "growth")) stop("argument should be a growth object")
  print(x@data$data)
  invisible(x)
}

