plotFreqB <- function(pos, freqB, pch=176, ylim=c(0,1), xlab="Position (Mb)", ylab=expression(beta == theta[B]/(theta[A]+theta[B]))) {
  xlim <- range(pos, na.rm=TRUE);
  scale <- 0.04*diff(xlim);
  xlim[1] <- xlim[1] - 3*scale;
  xlim[2] <- xlim[2] + 1*scale;
  x0 <- xlim[1];

  # Plot raw data
  plot(pos, freqB, pch=pch, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab);

  # Plot density
  d <- density(na.omit(freqB), from=0, to=1);
  d$y <- d$y / max(d$y, na.rm=TRUE);
  x <- d$x;
  d$x <- x0 + 2.5*scale*d$y;
  d$y <- x;
  lines(d, lwd=2, col="blue");

  # Annotate
#  abline(h=c(1/3,2/3), lty=2, lwd=2, col="blue");
  adj <- c(-0.1, 0.5);
  y <- 1/20;
  for (ss in c("AA", "AB", "BB")) {
    text(x=x0, y=y, adj=adj, ss, cex=1.3, col="blue");
    y <- y + 9/20;
  }
} # plotFreqB()


############################################################################
# HISTORY:
# 2009-01-27 [HB]
# o Created.
############################################################################
