plotdiversity <-
function(X, sample.times, makeplot=TRUE, filter=FALSE, ...) {
  if (filter) {
    seqn <- sample.times
  } else {
    seqn <- 1:length(sample.times)
  }
  diversity <- numeric(length(sample.times))
  k <- 1
  for (i in seqn) {
    diversity[k] <- meansnps(X$obs.strain[[i]], X$obs.freq[[i]], X$libr, X$nuc, X$librstrains)
    k <- k+1
  }
  if (makeplot) {
    plot(sample.times, diversity, type="l", ...)
  }
  return(invisible(diversity))
}
