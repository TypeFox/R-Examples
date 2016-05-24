setConstructorS3("CopyNumberOutliers", function(start=NULL, stop=NULL, mean=NULL, count=NULL, call=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  extend(Object(), "CopyNumberOutliers",
    start = start,
    stop = stop,
    mean = mean,
    count = count,
    call= call,
    ...
  )
})


setMethodS3("as.character", "CopyNumberOutliers", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, sprintf("Number of regions: %d", nbrOfRegions(this)));
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  GenericSummary(s);
}, protected=TRUE)


setMethodS3("nbrOfRegions", "CopyNumberOutliers", function(this, ...) {
  length(this$start);
})


setMethodS3("as.data.frame", "CopyNumberOutliers", function(x, ...) {
  # To please R CMD check
  this <- x;

  data <- cbind(start=this$start, stop=this$stop, mean=this$mean, count=this$count, call=this$call);
  data;
})


setMethodS3("applyRegions", "CopyNumberOutliers", function(this, FUN, ...) {
  data <- as.data.frame(this);
  o <- order(data[,"start"]);
  data <- data[o,,drop=FALSE];
  apply(data, MARGIN=1, FUN=FUN);
}, protected=TRUE)


setMethodS3("drawLevels", "CopyNumberOutliers", function(this, col="red", lwd=2, lty=1, xScale=1e-6, yScale=1, ...) {
  col0 <- col;
  lwd0 <- lwd;
  lty0 <- lty;
  applyRegions(this, FUN=function(cnr) {
    x <- c(cnr["start"], cnr["stop"]);
    y <- rep(cnr["mean"], times=2);
    if (is.function(col0))
      col <- col0(cnr);
    if (is.function(lwd0))
      lwd <- lwd0(cnr);
    if (is.function(lty0))
      lty <- lty0(cnr);
    lines(x=xScale*x, y=yScale*y, col=col, lwd=lwd, lty=lty, ...);
  })
})



setMethodS3("lines", "CopyNumberOutliers", function(x, col="red", lwd=2, xScale=1e-6, yScale=1, ...) {
  # To please R CMD check.
  this <- x;

  data <- as.data.frame(this);
  o <- order(data[,"start"]);
  data <- data[o,,drop=FALSE];
  xx <- t(data[,c("start", "stop"),drop=FALSE]);
  yy <- rep(this$mean[o], each=2);
  lines(x=xScale*xx, y=yScale*yy, col=col, lwd=lwd, ...);
})




############################################################################
# HISTORY:
# 2008-05-17
# o Moved to aroma.core.
# 2007-08-22
# o Created.  Need a generic container for holding copy number regions and
#   to plot them nicely.
############################################################################
