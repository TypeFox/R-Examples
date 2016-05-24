print.boostrWrapper <- function(x, ...) {
  attrs <- attributes(x)
  print(attrs[c("className", ".options", "predictor")])
}