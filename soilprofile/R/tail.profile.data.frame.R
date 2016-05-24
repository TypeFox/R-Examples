tail.profile.data.frame <-
function(x, ...) {
  class(x) <- 'data.frame'
  tail(x, ...)
}
