setMethodS3("extractC1C2", "list", function(fitList, ...) {
  c1c2List <- lapply(fitList, FUN=function(fit) {
    extractC1C2(fit, ...);
  });

  # Append NAs between chromosomes
  c1c2TList <- lapply(c1c2List, FUN=function(c1c2) {
    rbind(c1c2, NA);
  })

  c1c2 <- Reduce(rbind, c1c2TList);

  w <- sqrt(c1c2[,4]);
  w <- w / sum(w, na.rm=TRUE);
  w <- w / mean(w, na.rm=TRUE);
  c1c2 <- cbind(c1c2, w=w);

  class(c1c2) <- unique(c("C1C2", class(c1c2)));

  c1c2;
});

setMethodS3("plot", "C1C2", function(x, xlim=c(0,4), ylim=xlim, xlab=expression(C[1]), ylab=expression(C[2]), ...) {
  plot(NA, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim);
  abline(a=0, b=1, lty=3);
  grid();
  points(x, ...);
})

setMethodS3("points", "C1C2", function(x, cex=sqrt(x[,"w"])+1/8, ...) {
  x <- x[,c("C1","C2"),drop=FALSE];
  NextMethod("points", object=x, cex=cex);
})


setMethodS3("fitLoess2D", "C1C2", function(X, Y, ...) {
  fit <- fitLoessKD(X=X[,1:2,drop=FALSE], Y=Y[,1:2,drop=FALSE]);
  class(fit) <- c("Loess2DFit", class(fit));
  fit;
})

setMethodS3("normalizeLoess2D", "C1C2", function(X, ...) {
  XN <- X;
  XN[,1:2] <- normalizeLoessKD(X[,1:2], ...);
  XN;
}) # normalizeLoess2D()


##############################################################################
# HISTORY
# 2014-02-04
# o ROBUSTNESS: Now points() for C1C2 passes (modified) argument 'x' to
#   NextMethod() as 'object=x'.
# 2012-09-21 [HB]
# o Now extractDeltaC1C2() makes sure that there are splitters between
#   chromosomes and gaps.
# 2010-10-05 [HB]
# o Added extractPolarDeltaC1C2(), which returns [0,pi] angles and lengths.
# o Now extractDeltaC1C2() returns two types of change-point weights.
# 2010-09-19 [HB]
# o Created.
##############################################################################
