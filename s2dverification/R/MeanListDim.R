MeanListDim <- function(var, dims, narm = TRUE) {
  whichdims <- sort(dims, decreasing = TRUE)
  for (jdim in 1:length(whichdims)) {
    if (whichdims[jdim] <= length(dim(var))) {
      if (length(dim(var)) > 1) {
        var <- Mean1Dim(var, posdim = whichdims[jdim], narm = narm)
      } else {
        var <- mean(var, na.rm = narm)
      }
    }
  }
  # 
  #  Outputs
  # ~~~~~~~~~
  # 
  outvar <- var
}
