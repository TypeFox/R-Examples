multiKernDisplay <-
function (kern, spaceNum=0) {
  spacing = matrix("", spaceNum+1)
  cat(spacing);
  cat("Multiple output block kernel:\n")
  for(i in seq(along=kern$comp)) {
    cat(spacing)
    cat("Block ", i, "\n", sep="")
    kernDisplay(kern$comp[[i]], spaceNum)
  }
}
