setMethodS3("plotCoverageMap", "AromaUgpFile", function(this, chromosomes=getChromosomes(this), pch=15, cex=0.4, ..., yOffset=0, chrLabels=sprintf("Chr%02d", chromosomes), xlim=c(-10/xScale, max(this[,2L,drop=TRUE], na.rm=TRUE)), ylim=c(0, length(chromosomes)+1), xlab="Position", xScale=1e-6, add=FALSE) {
  chromosome <- this[,1L, drop=TRUE];
  x <- this[,2L, drop=TRUE];
  x <- xScale * x;

  yIdxs <- seq_along(chromosomes);
  

  if (!add) {
    plot(NA, xlim=xScale*xlim, ylim=ylim, xlab=xlab, ylab="", axes=FALSE);
    axis(side=1);
    usr <- par("usr");
    dx <- diff(usr[1:2]);
    x0 <- 0;
    yy <- yIdxs;
    xx <- rep(x0, times=length(yy));
    text(x=xx, y=yy, chrLabels, cex=0.8, pos=2, offset=1, xpd=TRUE);
  }
  
  for (cc in seq_along(chromosomes)) {
    chr <- chromosomes[cc];
    units <- (chromosome == chr);
    xx <- x[units];
    yy <- rep(yIdxs[cc]+yOffset, times=length(xx));
    points(xx, yy, pch=pch, cex=cex, ...);
  } # for (cc ...)
}) # plotCoverageMap()


############################################################################
# HISTORY:
# 2011-12-10
# o Added plotCoverageMap() for AromaUgpFile.
# o Created.
############################################################################
