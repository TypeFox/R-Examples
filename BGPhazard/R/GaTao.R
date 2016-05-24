GaTao <-
function(times, delta, type.t, K) {
  t.unc <- sort(times[delta  ==  1])
  if (type.t == 1) {
    n <- length(t.unc)
    if (n > K) {
      quant <- c(t.unc, max(times))
      tao <- c(0, quantile(x = times, probs = (1:K) / K, names = FALSE))
      if (type.t == 1 && length(unique(tao)) != length(tao)) {
        warning("Too many repeated observations. Zero-length intervals may
                appear.")
      }
    } else {
      stop (paste("The partition length (", K,") must be smaller than the number 
                  of uncensored times (", n, ").", sep = ""))
    }
  }
  if (type.t == 2) {
    K.t2 <- ceiling(max(times))
    if (K != K.t2) {
      aux <- K
      K <- K.t2       
      warning (c("'type.t' 2 requires K = ", K.t2, ". K (", aux,
                 ") fixed at ", K.t2, "."))
    }
    tao <- seq(0, K)
  }
  if (type.t  ==  3) {
    tao <- seq(0, ceiling(max(times)), ceiling(max(times)) / K)
  }
  return(tao)
}
