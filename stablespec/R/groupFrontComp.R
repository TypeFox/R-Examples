groupFrontComp <- function(matOfFronts, stringSize) {
  allComp <- list()
  b <- max(matOfFronts[, stringSize + 2])

  for (i in 0:b) {

    tempMatB <- matOfFronts[which(matOfFronts[, stringSize + 2] == i), 1:stringSize]

    if (length(tempMatB) > 0) {
      if(is.vector(tempMatB)) {
        tempMatB <- matrix(tempMatB, 1, stringSize, byrow=TRUE)
      }
      allComp[[i+1]] <- tempMatB
    } else {
      allComp[[i+1]] <- matrix(0, 1, stringSize)
    }
  }
  return(allComp)
}
