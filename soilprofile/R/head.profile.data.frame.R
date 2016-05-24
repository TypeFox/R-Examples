head.profile.data.frame <-
function(x, ...) {
  class(x) <- 'data.frame'
  head(x, ...)
}
