SKMtest <- function(y, trt, n, dferror, mserror, alpha, dms) {
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
  } else {
    Ymean <- mean(Ybar[(posmax + 1):n])
  }
  qobs <- (Ybar[1] + Ybar[n]) / 2 - Ymean
  groups <- rep(0, times = n)
  ng <- 1
  if ((qobs >= - dms[1]) & (qobs <= dms[1])) {
    groups[1:n] <- 1  
  } else {
    groups[1:posmax] <- 1
    groups[(posmax + 1):n] <- 2
    ng <- ng + 1
  }
  if ((any(groups != 1))) {
    continua1 = TRUE
  } else {
    continua1 = FALSE
  }
  if (continua1 == TRUE) {
    posI <- 1
    fim <- FALSE
    repeat {
      ini  <- groups[posI]
      posF <- max(which(groups == ini))
      if (((posF-posI) > 0) & (ini > 0)) {
        gap <- Ybar[(posI + 1):posF] - Ybar[posI:(posF - 1)]
        if(length(max(gap)==gap) > 1){
          valmax <- which(gap == max(gap))
          auxpos <- min(c(posF - valmax[1], valmax[1]))
          for(i in valmax[-1]){
            auxpos <- cbind(auxpos,min(c(posF - i, i)))
          }
          auxpos <- which.max(auxpos)
          posmax <- valmax[auxpos]
        } else{
          posmax <- as.integer(which.max(Ybar[(posI + 1):posF] 
                                         - Ybar[posI:(posF - 1)]))
        }
        if (posmax >= (posF - posmax)) {
          Ymean <- mean(Ybar[posI:(posI + posmax - 1)])
        } else {
          Ymean <- mean(Ybar[(posI + posmax):posF])
        }
        qobs <- (Ybar[posI] + Ybar[posF]) / 2 - Ymean
        aux <- groups[posI]
        if ((qobs >= - dms[2]) & (qobs <= dms[2])){
          groups[posI:posF] <- -aux
          posI <- posF + 1
          if (posI > n) {
            posI <- 1
          }
        } else {
          groups[(posI + posmax):posF] <- ng + 1
          ng <- ng + 1
          posI <- 1
        }
      } else{
        posI <- posF + 1
      }
      if (posI >= n) {
        fim <- TRUE
      }
      if (fim == TRUE) break
    }
    groups <- abs(groups)
  }
  result <- cbind(Ybar, groups)
  return(group.test2(result))
}