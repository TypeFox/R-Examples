#############################################################
#
#	plot.localdepth function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: February, 07, 2011
#	Version: 0.2
#
#	Copyright (C) 2011 Claudio Agostinelli and Mario Romanazzi
#
#############################################################

plot.localdepth <- function(x, xlab="Depth", ylab="Local Depth", main="DD plot", mark=0.9, labels=NULL, ...) {
  if (is.null(labels))
    labels <- rownames(x$y)
  sdepth <- x$depth
  ldepth <- x$localdepth
  if (max(sdepth)-min(sdepth))
    sdepth <- (sdepth-min(sdepth))/(max(sdepth)-min(sdepth))
  else
    sdepth <- rep(1, length(sdepth))
  if (max(ldepth)-min(ldepth))
    ldepth <- (ldepth-min(ldepth))/(max(ldepth)-min(ldepth))
  else
    ldepth <- rep(1, length(ldepth))
  plot(x=sdepth, y=ldepth, xlab=xlab, ylab=ylab, main=main, ...)
  tomark <- sdepth > mark | ldepth > mark
  text(sdepth[tomark], ldepth[tomark], labels=labels[tomark], pos=3)
  invisible(x)
}
