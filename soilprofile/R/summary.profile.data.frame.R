summary.profile.data.frame <-
function(object, ...) {
  class(object) <- 'data.frame'
  summary(object, ...)
}
