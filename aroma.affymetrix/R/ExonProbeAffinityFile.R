###########################################################################/**
# @RdocClass ExonProbeAffinityFile
#
# @title "The ExonProbeAffinityFile class"
#
# \description{
#  @classhierarchy
#
#  This class represents estimates of probe affinities in exon array
#  probe-level models.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "ProbeAffinityFile".}
#   \item{mergeGroups}{Specifies if the groups (exons) are merged or not for
#      these estimates.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "KS, HB"
#*/###########################################################################
setConstructorS3("ExonProbeAffinityFile", function(..., mergeGroups=FALSE) {
  extend(ProbeAffinityFile(...), "ExonProbeAffinityFile",
    mergeGroups=mergeGroups
  )
})


###########################################################################/**
# @RdocMethod getCellIndices
#
# @title "Retrieves tree list of cell indices for a set of units"
#
# \description{
#   @get "title" from the associated CDF.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Additional arguments passed to \code{getCellIndices()}
#             of @see "ProbeAffinityFile".}
#  \item{unlist}{If @TRUE, the cell indices are returned as a @vector.}
# }
#
# \value{
#   Returns a @list structure, where each element corresponds to a unit.
#   If argument \code{unlist=TRUE} is passed, an @integer @vector is returned.
# }
#
# \seealso{
#   @seeclass
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("getCellIndices", "ExonProbeAffinityFile", function(this, ..., unlist=FALSE) {
  # Argument 'unlist':
  unlist <- Arguments$getLogical(unlist);


  # Supported case?
  mergeGroups <- this$mergeGroups;
  if (unlist && mergeGroups) {
    throw("Unsupported request: Argument 'unlist' have to be TRUE when parameter 'mergeGroups' is TRUE: ", unlist);
  }

  cells <- NextMethod("getCellIndices");

  # Merge groups?
  if (mergeGroups) {
    cells <- .applyCdfGroups(cells, cdfMergeGroups);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Unlist?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (unlist) {
    cells <- unlist(cells, use.names=FALSE);
  }


  cells;
}, protected=TRUE) # getCellIndices()


setMethodS3("setMergeGroups", "ExonProbeAffinityFile", function(this, status, ...) {
  this$mergeGroups <- status;
}, protected=TRUE)


############################################################################
# HISTORY:
# 2007-02-07
# o Created (based on SnpProbeAffinityFile.R following chat with HB on
#   2007-02-07).
############################################################################
