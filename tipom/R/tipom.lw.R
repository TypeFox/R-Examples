tipom.lw <-
  function (lengths, widths, ia=FALSE, ms=FALSE, ...) {
    tipom.plot <- plot(widths, lengths,
                       ylim=c(0, max(lengths)+max(lengths)*0.1),
                       xlim=c(0, max(widths)+max(widths)*0.1),
                       asp=1,
                       xlab="Width",
                       ylab="Length",
                       ...)
    if (ms == TRUE) {
      abline(a=40, b=-1)
      abline(a=20, b=-1)
    }
    if (ia == TRUE) {
      slopes <- c(0.5, 0.75, 1, 1.5, 2, 3, 6)
      for (s in 1:7)
        abline(a=0, b=slopes[s])
    }
    tipom.plot
  }
