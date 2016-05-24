setMethodS3("drawC1C2Density", "PairedPSCBS", function(fit, ...) {
  # Nothing todo?
  if (nbrOfSegments(fit) < 2) {
    return(invisible());
  }

  data <- extractC1C2(fit);
  n <- data[,4, drop=TRUE];
  n <- sqrt(n);
  w <- n/sum(n, na.rm=TRUE);
  adjust <- 0.2;

  # For each dimension...
  for (cc in 1:2) {
    y <- data[,cc];
    ok <- is.finite(y) & is.finite(w);
    y <- y[ok];
    wt <- w[ok]/sum(w[ok]);

    # Nothing to do?
    if (length(y) < 2) {
      next;
    }

    d <- density(y, weights=wt, adjust=adjust);
    draw(d, side=cc, height=0.3, col="gray", lwd=2, xpd=FALSE);
    if (cc == 2) {
      draw(d, side=1, height=0.3, col="lightblue", lwd=2, xpd=FALSE);
    }
    p <- .findPeaksAndValleys(d, tol=0.05);
    type <- NULL; rm(list="type"); # To please R CMD check
    p <- subset(p, type == "peak");
    p <- p[order(p$density, decreasing=TRUE),,drop=FALSE];
    p <- head(p, n=8);
    if (cc == 1) {
      abline(v=p$x, lty=3, col="gray");
    } else {
      abline(h=p$x, lty=3, col="gray");
    }
  }
  box();
}) # drawC1C2Density()


setMethodS3("plotC1C2Grid", "PairedPSCBS", function(fit, ..., Clim=c(0,4), main=NULL) {
  plotC1C2(fit, ..., Clim=Clim);
  title(main=main);
  drawC1C2Density(fit, ...);
})



##############################################################################
# HISTORY
# 2012-04-16
# o CLEANUP: drawC1C2Density() for PairedPSCBS no longer needs to require
#   'aroma.core', because draw() for 'density objects are in R.utils 1.10.0.
# 2012-02-27
# o BUG FIX: drawC1C2Density() for PairedPSCBS would throw an exception
#   if there was only one segment, or less than two finite (C1,C2):s.
# 2012-02-23
# o Moved drawC1C2Density() and plotC1C2Grid() to its own source file.
# 2011-11-12 [HB]
# o Added drawC1C2Density() adopted from plotC1C2Grid().
# o Now argument '...' of plotC1C2Grid() for PairedPSCBS are passed
#   to plotC1C2().
# 2011-07-10 [HB]
# o Updated code to work with the new column names in PSCBS v0.11.0.
# 2010-10-08 [HB]
# o Added plotC1C2Grid() for PairedPSCBS.
# o Created.
##############################################################################
