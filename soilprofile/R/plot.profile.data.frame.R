plot.profile.data.frame <-
function(x, y, ...) {
  class(x) <- 'data.frame'
  plot_profile(x, ...)
}
