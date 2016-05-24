###########################################################################/**
# @set "class=AbstractCBS"
# @RdocMethod plotTracks
#
# @title "Plots the segmentation result along the genome"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{...}
# }
#
# \value{
#   Returns nothing.
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("plotTracks", "AbstractCBS", abstract=TRUE);


setMethodS3("tileChromosomes", "AbstractCBS", abstract=TRUE, protected=TRUE);


setMethodS3("drawChangePoints", "AbstractCBS", abstract=TRUE, protected=TRUE);


setMethodS3("drawKnownSegments", "AbstractCBS", function(fit, col="#aaaaaa", ..., xScale=1e-6) {
  segs <- fit$params$knownSegments;

  # Nothing todo?
  if (is.null(segs)) {
    return();
  }

  # Workaround from the fact that extractChromosomes() does not drop
  # known segments. /HB 2013-03-21
  chromosome <- NULL; rm(list="chromosome");  # To please R CMD check.
  segs <- subset(segs, chromosome %in% getChromosomes(fit));
  xStarts <- segs[,"start"];
  xEnds <- segs[,"end"];
  xs <- sort(unique(c(xStarts, xEnds)));
  abline(v=xScale*xs, col=col, ...);
}, protected=TRUE)


############################################################################
# HISTORY:
# 2013-03-21
# o Added drawKnownSegments().
# 2011-12-03
# o Added drawChangePoints().
# 2011-10-02
# o Created.
############################################################################
