poslin <-
function(x) {
  y <- x
  y[y < 0] <- 0
  y
}
