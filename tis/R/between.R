between <- function(y, x1, x2){
  y <- unclass(y)
  x1 <- unclass(x1)
  x2 <- unclass(x2)
  small <- pmin(x1, x2)
  large <- pmax(x1, x2)
  (y >= small) & (y <= large)
}
