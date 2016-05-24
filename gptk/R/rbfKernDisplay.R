rbfKernDisplay <-
function (kern, spaceNum=0) {
  spacing <- matrix("", spaceNum+1)
  cat(spacing)
  if ("isNormalised" %in% names(kern) && kern$isNormalised)
    cat("Normalised version of the kernel.\n")
  else
    cat("Unnormalised version of the kernel.\n")
  cat(spacing)
  cat("RBF inverse width: ", kern$inverseWidth, " (length scale ", 1/sqrt(kern$inverseWidth), ")\n", sep="")
  cat(spacing)
  cat("RBF variance: ", kern$variance, "\n", sep="")
}
