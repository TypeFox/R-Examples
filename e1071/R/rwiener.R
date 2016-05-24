rwiener <- function(end=1, frequency=1000) {

  z<-cumsum(rnorm(end*frequency)/sqrt(frequency))
  ts(z, start=1/frequency, frequency=frequency)
}

