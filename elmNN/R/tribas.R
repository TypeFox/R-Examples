tribas <-
function(x) {
  y <- ifelse(x >= -1 & x <= 1, 1 - abs(x), 0)
  y
}
