mid.in <- function(x,y,rx,ry) {
xm <- 0.5*(x[1]+x[2])
ym <- 0.5*(y[1]+y[2])
(rx[1] < xm & xm < rx[2] & ry[1] < ym & ym < ry[2])
}
