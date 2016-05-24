###########################################################################/**
# @set "class=list"
# @RdocMethod plotBoxplotStats
#
# @title "Plots a (merged or non-merged) list of boxplot.stats() elements"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{stats}{A (merged or non-merged) @list of
#                                  @see "grDevices::boxplot.stats" elements.}
#   \item{ylim, outline, las, ...}{Arguments passed to @see "graphics::bxp".}
# }
#
# \value{
#   Returns (invisibly) the merged @see "grDevices::boxplot.stats" structure.
# }
#
# @author "EP"
#
# @keyword internal
#*/###########################################################################
setMethodS3("plotBoxplotStats", "list", function(stats, ylim=NULL, outline=FALSE, las=2, ...) {
  bxpStats <- mergeBoxplotStats(stats);

  # fix the strange behavior of bxp if outline=FALSE
  if(is.null(ylim) && !outline) {
    ylim <- range(as.vector(bxpStats[["stats"]]));
  }

  bxp(bxpStats, ylim=ylim, outline=outline, las=las, ...);

  invisible(stats);
})

##########################################################################
# HISTORY:
# 2008-02-25
# o Renamed.
# 2008-02-22 [EP]
# o Created.
##########################################################################

