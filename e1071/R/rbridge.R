rbridge <- function(end=1, frequency=1000) {

  z <- rwiener(end=end, frequency=frequency)
  ts(z - time(z)*as.vector(z)[frequency],
     start=1/frequency, frequency=frequency)
}

