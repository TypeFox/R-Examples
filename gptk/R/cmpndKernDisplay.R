cmpndKernDisplay <-
function (kern, spaceNum=0) {
  spacing = matrix("", spaceNum+1)
  cat(spacing)
  cat("Compound kernel:\n")
  for(i in seq(along=kern$comp)) 
    kernDisplay(kern$comp[[i]], spaceNum+2)
}
