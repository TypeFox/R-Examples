`head.growth` <-
function(x, n = 6, ...) {
  if (!is(x, "growth")) stop("argument should be a growth object")
  xh = x@data$data[1:min(n, nrow(x@data$data)),]
  print(xh)
  invisible(x)
}

