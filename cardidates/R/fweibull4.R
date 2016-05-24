`fweibull4` <-
function(x, p) {
  ## dweibull does not like named vectors as first argument
  p[1] + p[2] * dweibull(as.numeric(x), p[3], p[4])
}

