a_get_relative <-
function(x, y){
coords <-par('usr')
dinps  <-par('pin')
xr <- abs(diff(coords[1:2]))
yr <- abs(diff(coords[3:4]))
xtox   <- dinps[1]/dinps[2]
x <- (x-coords[1])/xr
y <- (y-coords[3])/yr
x  <- x*xtox
cbind(x, y)
}
