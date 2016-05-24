build.hull.list <-
function(frame) {
  ## Given a data frame containing isotope values,
  ## return a matrix containing the values in the data frame, "stacked" for
  ## the use of convhulln
  d = dim(frame)[2]
  out = c()
  for (i in seq(1, length(names(frame)), d)) {
    a = frame[i:(i+d-1)]
    a = na.omit(a)
    out = rbind(out, a)
  }
  out
}
