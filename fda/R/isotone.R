isotone <- function(y) {
	#  Compute an isotonic regression line from data in Y
	#  This is the piecewise linear line that is monotonic and 
	#  most closely approximates Y in the least squares sense.
  n    <- length(y)
  mony <- y
  eb   <- 0
  indx <- 1:(n-1)
  while (eb < n) {
    negind <- (diff(mony) < 0)
    ib <- min(indx[diff(mony) < 0])
    if (is.na(ib)) {
      bb <- eb <- n
    } else {
      bb <- eb <- ib
    } 
    while (eb < n && mony[bb] == mony[eb+1]) eb <- eb + 1
    poolflg <- -1
    while (poolflg != 0) {
      if (eb >=  n || mony[eb] <= mony[eb+1]) poolflg <- 1
      if (poolflg == -1) {
        br <- er <- eb+1
        while (er < n && mony[er+1] == mony[br]) er <- er + 1
        pmn <- (mony[bb]*(eb-bb+1) + mony[br]*(er-br+1))/(er-bb+1)
        eb <- er
        mony[bb:eb] <- pmn
        poolflg <- 1
      }
      if (poolflg == 1) {
        if (bb <= 1 || mony[bb-1] <= mony[bb]) {
          poolflg <- 0
        } else {
          bl <- el <- bb-1
          while (bl > 1 && mony[bl-1] == mony[el]) bl <- bl - 1
          pmn <- (mony[bb]*(eb-bb+1) + mony[bl]*(el-bl+1))/(eb-bl+1)
          bb <- bl
          mony[bb:eb] <- pmn
          poolflg <- -1
        }
      }
    }
  }
  return(mony)
}
