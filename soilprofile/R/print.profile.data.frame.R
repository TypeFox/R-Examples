print.profile.data.frame <-
function(x, ...) {
  class(x) <- 'data.frame'
  print(x, ...)
}
