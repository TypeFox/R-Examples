bijectivityback <-
function(xt, power)
 {
  if (power == 0)
   {
   	res <- T } else {
   	x <- powerUntransform(xt = xt, power = power)
    newxt <- powerTransform(x = x, power = power)
    res <- all.equal(xt, newxt)
   }
 }
