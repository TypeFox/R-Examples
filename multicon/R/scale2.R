scale2 <-
function(x, center = TRUE, scale = TRUE) {
  x <- data.frame(x)
  miss <- colSums(is.na(x))
  valid <- nrow(x) - miss
  x.means <- colMeans(x, na.rm=T)
  x.sds <- sqrt((apply(x, 2, sd, na.rm=T)^2) * (valid-1) / (valid))

  if(center==T & scale==T) {
    out <- t((t(x) - x.means) / x.sds)
  }
  if(center==T & scale==F) {
    out <- t(t(x) - x.means)
  }
  if(center==F & scale==T) {
    out <- t(t(x) / x.sds)
    warning("Standardizing without centering is unconventional.")
  }
  if(center==F & scale==F) {
    out <- x
    warning("You realize you didn't do anything to your data right?")
  }
  return(out)
}
