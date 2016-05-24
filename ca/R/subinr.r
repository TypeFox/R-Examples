################################################################################
# subinr(): Computing inertia of 'sub-matrices' (ca package 0.70)
################################################################################
subinr <- function(B, ind) {
  nn   <- length(ind)
  subi <- matrix(NA, nrow = nn, ncol = nn)
  ind2 <- c(0,cumsum(ind))
  for (i in 1:nn) {
    for (j in 1:nn) {
      tempmat   <- B[(ind2[i]+1):(ind2[i+1]), (ind2[j]+1):(ind2[j+1])]
      tempmat   <- tempmat / sum(tempmat)
      er        <- apply(tempmat, 1, sum)
      ec        <- apply(tempmat, 2, sum)
      ex        <- er%*%t(ec)
      subi[i,j] <- sum((tempmat - ex)^2 / ex)
      }
    }
  return(subi/nn^2)
  }
################################################################################
