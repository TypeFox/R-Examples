setMethodS3("drawLevels", "DNAcopy", function(fit, field=c("seg.mean", "tcn.mean", "dh.mean"), xScale=1, col="red", lwd=3, ...) {
  field <- match.arg(field);
  segments <- fit$output[,c("loc.start", "loc.end", field)];
  apply(segments, MARGIN=1, FUN=function(seg) {
    x <- c(seg[["loc.start"]], seg[["loc.end"]]);
    y <- rep(seg[[field]], times=2);
    lines(x=xScale*x, y=y, col=col, lwd=lwd, ...);
  });
})



############################################################################
# HISTORY:
# 2010-07-09
# o Created from drawLevels() for CopyNumberRegions in aroma.core.
############################################################################
