plotccda.results <-
function (x) 
  {
    ymin=0
    if(min(x$difference)<0){ymin=min(x$difference)*110}
    ngroups = length(x$ratio)
    par(mfrow = c(1, 2))
    plot(x$ratio * 100, type = "o", col = "black", ylim = c(ymin, 
                                                            130), lwd = 2, xlab = "number of groups", ylab = "LDA-percentages (%)", 
         main = "Results", font.main = 3
         , xaxt = "n"
    )
    axis(1, at = 1:ngroups)
    lines(x$q95 * 100, type = "o", col = "green", lty = 2, lwd = 2)
    lines(x$difference * 100, type = "o", col = "blue", lty = 3, 
          lwd = 2)
    legend("top", c("ratio", "q95", "d=ratio-q95"), cex =min(0.9*par('pin')[1]/strwidth("d=ratio-q95","inches"),par('pin')[2]*(0.4/13)/strheight("d=ratio-q95","inches")), 
           col = c("black", "green", "blue"), pch = c(1, 1, 1), 
           lty = c(1, 2, 3), lwd = c(2, 2, 2))
    plot(x$difference * 100, type = "o", main = "Difference", 
         col = "blue", lty = 3, lwd = 2, xlab = "number of groups", 
         ylab = "LDA-differences (%)", font.main = 3, xaxt = "n")
    axis(1, at = 1:ngroups)
  }
