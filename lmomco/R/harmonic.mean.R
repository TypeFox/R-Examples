"harmonic.mean" <-
function(x) {
  x <- x[! is.na(x)]
  n <- length(x)
  x.nonzero <- x[x != 0]
  n.zero <- n - length(x.nonzero)
  correction <- (n-n.zero)/n
  HM <- (1/mean(1/x.nonzero)) * correction
  z <- list(harmean = HM, 
            correction = correction,
            source = "harmonic.mean")
  return(z)
}
