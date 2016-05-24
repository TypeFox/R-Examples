###########################################################################/**
# @set "class=list"
# @RdocMethod mergeBoxplotStats
#
# @title "Merges a list of boxplot.stats() elements"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{stats}{A @list of elements each in a format returned by
#      @see "grDevices::boxplot.stats".}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns a @list structure in a format that is returned by
#   @see "graphics::boxplot".
# }
#
# @examples "../incl/mergeBoxplotStats.list.Rex"
#
# @author
#
# @keyword internal
#*/###########################################################################
setMethodS3("mergeBoxplotStats", "list", function(stats, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  isBoxplotStats <- function(bxp, ...) {
    if (!is.list(bxp))
      return(FALSE);

    if (!all(c("stats", "conf") %in% names(bxp)))
      return(FALSE);

    TRUE;
  } # isBoxplotStats()



  # Do nothing?
  if (isBoxplotStats(stats))
    return(stats);


  # Sanity check
  for (kk in seq_along(stats)) {
    stat <- stats[[kk]];
    if (!isBoxplotStats(stat)) {
      throw("Element #", kk, " in argument 'stats' is not a list structure returned by boxplot.stats(): ", class(stat)[1]);
    }
  }


  # Append 'group' stats
  for (kk in seq_along(stats)) {
    stats[[kk]][["group"]] <- rep(kk, times=length(stats[[kk]][["out"]]));
  }


  # Restructure so it is recognized by graphics::bxp().
  bxpStats <- list();

  for (field in names(stats[[1]])) {
    args <- lapply(stats, FUN=.subset2, field);
    value <- do.call(c, args);
    if (field == "stats") {
      value <- matrix(value, nrow=5);
    } else if (field == "conf") {
      value <- matrix(value, nrow=2);
    }
    names(value) <- NULL;
    dimnames(value) <- NULL;
    bxpStats[[field]] <- value;
  }

  bxpStats[["names"]] <- names(stats);

  bxpStats;
}, protected=TRUE)


##########################################################################
# HISTORY:
# 2008-02-25
# o Now it returned the input object if already in the correct format.
# o Renamed. Added a sanity check. Add Rdoc comments.
# 2008-02-22 [EP]
# o Created.
##########################################################################
