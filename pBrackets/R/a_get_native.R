a_get_native <-
function(x, y){
coords <-par('usr')
dinps  <-par('pin')
xtox   <- dinps[1]/dinps[2]
xr <- abs(diff(coords[1:2]))
yr <- abs(diff(coords[3:4]))
x  <- x/xtox
x <- ((x*xr)+coords[1])
y <- ((y*yr)+coords[3])
cbind(x, y)
}
