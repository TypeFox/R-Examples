setMethodS3("pairs2", "principal.curve", function(fit, pch=19L, cex=0.8, fitCol="red", fitLwd=2, fitLty=1L, xlim=NULL, ylim=xlim, lower.panel=NULL, ...) {
  r <- range(c(fit$s, fit$Y), na.rm=TRUE);

  # Argument 'xlim' & 'ylim':
  if (is.null(xlim)) {
    xlim <- r;
  }
  if (is.null(ylim)) {
    ylim <- r;
  }

  hasData <- !is.null(fit$Y);

  nbrOfVars <- ncol(fit$s);
  layout(matrix(seq_len((nbrOfVars-1)^2), nrow=nbrOfVars-1L, ncol=nbrOfVars-1L, byrow=TRUE));
  par(mar=c(3,3,1,1)+0.1);
  for (rr in seq_len(nbrOfVars)) {
    if (rr == nbrOfVars)
      next;
    for (cc in seq_len(nbrOfVars)) {
      if (cc < rr) {
        plot.new();
      } else if (cc == rr) {
      } else {
        plot(NA, xlim=xlim, ylim=ylim, xlab="", ylab="");
        abline(a=0, b=1, lty=3L, col="#999999", lwd=2);
        if (hasData) {
          y <- fit$Y[,c(cc,rr),drop=FALSE];
          points(y, pch=pch, cex=cex, ...);
        }

        if (fitLwd > 0) {
          y <- fit$s[,c(cc,rr),drop=FALSE];
          o <- order(y[,1L,drop=TRUE]);
          y <- y[o,,drop=FALSE];
          lines(y, col=fitCol, lwd=fitLwd, lty=fitLty);
        }
      }
    } # for (cc ...)
  } # for (rr ...)
}) # pairs2()


###########################################################################
# HISTORY:
# 2013-10-08
# o BUG FIX: pairs2() for 'principal.curve' would assume that the
#   fitted curve was ordered by the first dimension in each panel.
# 2009-01-12
# o Extracted from MultiSourceCopyNumberNormalization.R.
###########################################################################
