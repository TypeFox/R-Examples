powerUntransform <-
function(xt, power)
 {
  if (power == 0) x <- exp(xt) else x <- xt^(1/power)
  x
 }
