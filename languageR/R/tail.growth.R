`tail.growth` <-
function(x, n = 6, ...) {
  if (!is(x, "growth")) stop("argument should be a growth object")
  m = nrow(x@data$data)
  xh = x@data$data[max((m-n+1),0):m,]
  print(xh)
  invisible(x)
}

