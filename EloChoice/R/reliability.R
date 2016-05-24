# reliability 15_10_07
# replaces 'upsetindex()' function

reliability <- function(x) {
  ups <- N <- upsw <- numeric(as.numeric(x$misc["runs"]))
  i = 1
  for (i in 1:length(ups)) {
    ind <- which(x$decmat[i, ])
    u <- as.numeric(x$upsmat[i, ][ind])
    wg <- x$wgtmat[i, ][ind]
    ups[i] <- 1 - sum(u)/length(u)
    upsw[i] <- 1 - weighted.mean(u, w = abs(wg))
    N[i] <- length(u)
  }
  return(data.frame(upset = ups, upset.wgt = upsw, totIA = N))
}



