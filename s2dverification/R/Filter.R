Filter <- function(xdata, freq) {
  fac1 <- 1
  fac2 <- 1
  ndat <- length(xdata)
  ndat2 <- length(xdata[is.na(xdata) == FALSE])
  maxi <- 0
  endphase <- 0
  for (jfreq in seq(freq - 0.5 / ndat2, freq + 0.5 / ndat2, 0.1 / (ndat2 * 
       fac1))) {
    for (phase in seq(0, pi, (pi / (10 * fac2)))) {
      xtest <- cos(phase + c(1:ndat) * jfreq * 2 * pi)
      test <- lm(xdata[is.na(xdata) == FALSE] ~ xtest[
              is.na(xdata) == FALSE])$fitted.value
      if (sum(test ^ 2) > maxi) { 
        endphase <- phase
        endfreq <- jfreq
      }
      maxi <- max(sum(test ^ 2), maxi)
    }
  }
  xend <- cos(endphase + c(1:ndat) * endfreq * 2 * pi)
  xdata[is.na(xdata) == FALSE] <- xdata[is.na(xdata) == FALSE] - lm(
                              xdata[is.na(xdata) == FALSE] ~ xend[is.na(xdata) == FALSE]
                              )$fitted.value

  xdata
}
