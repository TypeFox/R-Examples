#### Maximum draw down ####

DrawDown <- function (x) {
  x <- na.omit(x)
  max.dd = 0
  if (length(x) < 2) {
    return (NA)
  }
  for (i in 1:length(x)) {
    temp.dd = 1 - x[i] / max(x[1:i])
    if (temp.dd > max.dd) {
      max.dd = temp.dd
    }
  }
  return (max.dd)
}