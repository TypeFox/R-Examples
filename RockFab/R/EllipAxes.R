EllipAxes <-
function(es, nu){
  y <- exp(1)^((es * nu * sqrt(6)) / (3 * sqrt(nu^2 + 3)))
  if(y != 1){
    x <- exp(1)^((log(y) * (3 - nu)) / (2 * nu))
    z <- exp(1)^((log(y) * (3 + nu)) / (-2 * nu))
  } else{
    x <- exp(1)^((es) / (sqrt(2)))
z <- 1 / x
  }
  return(c(x, y, z))
}
