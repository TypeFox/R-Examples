SelIndices <- function(var, posdim, limits) {
  #
  #  A few security checks
  # ~~~~~~~~~~~~~~~~~~~~~~~
  #
  dimsvar <- dim(var)
  if (is.null(dimsvar)) {
    dimsvar <- length(var)
  }
  if (posdim > length(dimsvar)) {
    stop("posdim does not exist")
  }
  if (length(limits) != 2) {
    stop("Need lower and upper limit")
  } 
  if (dimsvar[posdim] < limits[2]) {
    stop("Check the consistency between limits and var dimensions")
  }
  #
  #  Select the correct indices 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  enlvar <- Enlarge(var, 10)
  u <- IniListDims(dimsvar, 10)
  u[[posdim]] <- limits[1]:limits[2]
  enlvarout <- enlvar[u[[1]], u[[2]], u[[3]], u[[4]], u[[5]], u[[6]], u[[7]],
                      u[[8]], u[[9]], u[[10]]]
  #
  #  Preparing the output matrice
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  if (limits[2] == limits[1]) {
    dimsvar <- dimsvar[-posdim]
  } else {
    dimsvar[posdim] <- limits[2] - limits[1] + 1
  }
  varout <- array(dim = dimsvar)
  varout[] <- enlvarout
  
  #
  #  Output
  # ~~~~~~~~
  #
  varout
}
