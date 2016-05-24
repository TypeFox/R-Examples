oneByoneHack <- function(info)
  {
    info$rowsep <- numeric(0)
    info$colsep <- numeric(0)
    info$nrows <- 1
    info$ncols <- 1
    info$lengths <- info$lengths[1]
    info$rowReps <- info$rowReps[1]
    info$locations <- matrix(info$locations[1,],1,)
    info$indices <- info$indices[1:3,]
    info$indices.preLimits$ax <- info$indices.preLimits$ax[1]
    info$indices.preLimits$bx <- info$indices.preLimits$bx[1]
    info$indices.preLimits$getV.i <- info$indices.preLimits$getV.i[1]
    info$indices.preLimits$getV.j <- info$indices.preLimits$getV.j[1]
    info$indices.preLimits$indices <- matrix(info$indices.preLimits$indices[1,],1,)
    info$indices.preLimits$evalFactor <- info$indices.preLimits$evalFactor[1]
    info
  }
