diagReplace <- function(R, upper=TRUE) {
 RT <- R
 if (upper == TRUE) {
  Rtranspose         <- t(RT)
  # Replacing upper diagonal with lower diagonal
  RT[upper.tri(RT)]  <- Rtranspose[upper.tri(Rtranspose)]
  return(RT)
  }
 if (upper == FALSE) {
  Rtranspose         <- t(RT)
  # Replacing lower diagonal with upper diagonal
  RT[lower.tri(RT)]  <- Rtranspose[lower.tri(Rtranspose)]
  return(RT)
  }
 }


