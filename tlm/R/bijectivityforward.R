bijectivityforward <-
function(x, power)
 {
  if (power == 0)
   {
   	res <- !any(x <= 0) } else {
    xt <- powerTransform(x = x, power = power)
    newx <- powerUntransform(xt = xt, power = power)
    res <- all.equal(xt, newx)
   }
 }
