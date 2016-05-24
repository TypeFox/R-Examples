powerTransform <-
function(x, power)
 {
  if (power == 0) xt <- log(x) else xt <- x^power
  xt
 }
