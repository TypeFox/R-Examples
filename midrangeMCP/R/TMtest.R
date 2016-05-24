TMtest <- function(y, trt, n, dferror, mserror, alpha, dms) {
  Ybar <- tapply(y, trt, mean)
  Ybar <- sort(Ybar)
  gap <- Ybar[2:n] - Ybar[1:(n - 1)]
  if(length(max(gap)==gap) > 1){
    valmax <- which(gap == max(gap))
    auxpos <- min(c(n - valmax[1], valmax[1]))
    for(i in valmax[-1]){
      auxpos <- cbind(auxpos,min(c(n - i, i)))
    }
    auxpos <- which.max(auxpos)
    posmax <- valmax[auxpos]
  } else{
    posmax <- as.integer(which.max(Ybar[2:n] - Ybar[1:(n - 1)]))
  }
  if (posmax >= (n - posmax)) {
    Ymean <- mean(Ybar[1:posmax])
  } else
  {
    Ymean <- mean(Ybar[(posmax + 1):n])
  }
  range <- n
  pos <- 1
  col <- 0
  qobs <- (Ybar[pos] + Ybar[n]) / 2 - Ymean
  aux <- rep(0, times = n)
  groups <- Ybar
  if ((qobs >= - dms[1]) & (qobs <= dms[1])) aux[1:n] <- 1
  if (!(any(aux == 0))) {
    groups <- cbind(Ybar, aux)
  }
  if (any(aux == 0)) {
    continua1 = TRUE
  } else {
    continua1 = FALSE
  }
  if (continua1 == TRUE) {
    range <- range - 1
    pos <- 0
    ncomp <- n - range + 1
    ct <- 1
    if (range < 2) {
      continua2 <- FALSE
    } else {
      continua2 <- TRUE
    }
    if (continua2 == TRUE) {
      repeat {
        pos <- pos + 1
        posmax <- as.integer(which.max(Ybar[(pos + 1):(pos + range - 1)]
                                       - Ybar[pos:(pos + range - 2)]))
        gap <- Ybar[(pos + 1):(pos + range - 1)] - Ybar[pos:(pos + range - 2)]
        if(length(max(gap)==gap) > 1){
          valmax <- which(gap == max(gap))
          auxpos <- min(c(range - valmax[1], valmax[1]))
          for(i in valmax[-1]){
            auxpos <- cbind(auxpos,min(c(n - i, i)))
          }
          auxpos <- which.max(auxpos)
          posmax <- valmax[auxpos]
        }else{
          posmax <- as.integer(which.max(Ybar[(pos + 1):(pos + range - 1)]
                                         - Ybar[pos:(pos + range - 2)]))
        }
        if (posmax >= (pos + range - 1 - posmax)) {
          Ymean <- mean(Ybar[pos:(pos+posmax-1)])
        } else {
          Ymean <- mean(Ybar[(pos + posmax):(pos + range - 1)])
        }
        qobs <- (Ybar[pos] + Ybar[pos+range-1]) / 2 - Ymean
        aux[1:n] <- 0
        if (((qobs >= - dms[2]) & (qobs <= dms[2]))) {
          aux[pos:(pos + range - 1)] <- 1
        }
        dentro <- FALSE
        if ((col > 0) & any(aux == 1)) {
          for (i in 1:col)
          {
            if (any(groups[aux == 1, i + 1] == 0)) {
              dentro <- FALSE
            } else {
              dentro <- TRUE
            }
            if (dentro == TRUE) break
          }
        }
        if ((dentro == FALSE) & any(aux == 1)) {
          groups <- cbind(groups, aux)
          col <- col + 1
        }
        ct <- ct + 1
        if (ct > ncomp)
        {
          range <- range - 1
          pos <- 0
          ncomp <- n - range + 1
          ct <- 1
        }
        if (range < 2) {
          continua2 <- FALSE
        } else {
          continua2 <- TRUE
        }
        if (continua2 == FALSE) break
      }
    }
  }
  result <- ProcTest(groups)
  return(group.test(result))
}
