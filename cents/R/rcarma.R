#Source: rcarma.R
#function for simulating ARMA left or right censored time series
#allows multiple censoring points
#
rcarma <- function(n=200, ar=0.9, ma=0.6, mu=100, siga=15, rates = c(0.15, NA), Mrate=0) {
#rates - left, right with NA=not applicable
  Rates <- rates
  if (is.vector(rates))
    Rates <- matrix(rep(rates, n), byrow=TRUE, ncol=2)
  y <- z <- mu+siga*as.vector(arima.sim(model=list(ar=ar, ma=ma), n=n))
  iy <- yL <- yR <-rep(NA, n)
#
#left-censored
  cL <- quantile(z, Rates[,1])
  indL0 <- z > cL #not left-censored
  indL <- !ifelse(is.na(indL0), TRUE, indL0) #is left-censored
  y <- ifelse(indL, cL, z)  
#
#right-censored
  cR <- quantile(z, 1-Rates[,2])
  indR0 <- z < cR #not right-censored
  indR <- !ifelse(is.na(indR0), TRUE, indR0) #is right-censored
  y <- ifelse(indR, cR, y)  
#
#missing values
  indMissing <- is.element(1:n, sample(1:n, size=floor(Mrate*n)))
#some previously censored values may become missing!
  y[indMissing] <- yL[indMissing] <- yR[indMissing] <- NA
  indL <- indL & !indMissing
  indR <- indR & !indMissing
  indo <- !(indMissing|indL|indR)
#iy=o,L,R,NA according as fully observed, left/right censored, missing
  iy <- rep("na",n)
  iy[indo] <- "o"
  iy[indL] <- "L"
  iy[indR] <- "R"
#
  ans <- list(y=y, iy=iy, censorPts=matrix(c(cL,cR), ncol=2), z=z)
  class(ans) <- "cents"
  ans
}


