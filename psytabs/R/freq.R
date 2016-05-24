freq <-
  function (x) {
    freq <- data.frame(n = summary(x), pc = round(summary(x)/length(x)*100, 1))
    freq <- rbind(rep(NA, ncol(freq)), freq)
    row.names(freq) <- c(deparse(substitute(x)), paste("  ", levels(x), sep=""))
    freq
  }