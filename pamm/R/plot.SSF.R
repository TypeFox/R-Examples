`plot.SSF` <-
  function (x, ...) 
  {
    if (!inherits(x, "SSF")) 
      stop("use only with \"SSF\" xs")
    par(mfrow = c(2, 2))
    plot(x$nb.ID, x$int.pval, type = "l", xlab = "Nb Gr./Repl.", 
         ylab = "P-value", main = "Intercept P-value", ylim = c(0, 
                                                                1), xaxt = "n")
    axis(1, at = x$nb.ID, labels = paste(x$nb.ID, x$nb.repl, 
                                         sep = "/"), tick = FALSE)
    abline(h = 0.05)
    lines(x$nb.ID, x$CIup.ipv, lty = 2)
    lines(x$nb.ID, x$CIlow.ipv, lty = 2)
    plot(x$nb.ID, x$sl.pval, type = "l", xlab = "Nb Gr./Repl.", 
         ylab = "P-value", main = "Slope P-value", ylim = c(0, 
                                                            1), xaxt = "n")
    axis(1, at = x$nb.ID, labels = paste(x$nb.ID, x$nb.repl, 
                                         sep = "/"), tick = FALSE)
    abline(h = 0.05)
    lines(x$nb.ID, x$CIup.slpv, lty = 2)
    lines(x$nb.ID, x$CIlow.slpv, lty = 2)
    plot(x$nb.ID, x$int.power, type = "l", ylim = c(0, 1), xlab = "Nb Gr./Repl.", 
         ylab = "Power", main = "Intercept Power Calculations", 
         xaxt = "n")
    axis(1, at = x$nb.ID, labels = paste(x$nb.ID, x$nb.repl, 
                                         sep = "/"), tick = FALSE)
    lines(x$nb.ID, x$CIup.ipo, lty = 2)
    lines(x$nb.ID, x$CIlow.ipo, lty = 2)
    plot(x$nb.ID, x$sl.power, type = "l", ylim = c(0, 1), xlab = "Nb Gr./Repl.", 
         ylab = "Power", main = "Slope Power Calculations", xaxt = "n")
    axis(1, at = x$nb.ID, labels = paste(x$nb.ID, x$nb.repl, 
                                         sep = "/"), tick = FALSE)
    lines(x$nb.ID, x$CIup.slpo, lty = 2)
    lines(x$nb.ID, x$CIlow.slpo, lty = 2)
  }
