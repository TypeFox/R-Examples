satlins <-
function(x) {
  y <- ifelse(x >= 1, 1, ifelse(x <= -1, -1, x))
  y
}
