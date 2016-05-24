cmyk2rgb <- function(cmyk){
  R <- (cmyk[, 1L] - 1.0) * (cmyk[, 4L] - 1.0)
  G <- (cmyk[, 2L] - 1.0) * (cmyk[, 4L] - 1.0)
  B <- (cmyk[, 3L] - 1.0) * (cmyk[, 4L] - 1.0)

  return(rgb(R,G,B))
}
