whiteKernDisplay <-
function (kern, spaceNum=0) {
  spacing <- matrix("", spaceNum+1)
  cat(spacing)
  cat("White Noise Variance: ", kern$variance, "\n", sep="")
}
