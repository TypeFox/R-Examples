###########################################################################/**
# @RdocClass SnpProbeAffinityFile
#
# @title "The SnpProbeAffinityFile class"
#
# \description{
#  @classhierarchy
#
#  This class represents estimates of probe affinities in SNP probe-level models.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "ProbeAffinityFile".}
#   \item{mergeStrands}{Specifies if the strands are merged or not for these
#      estimates.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("SnpProbeAffinityFile", function(..., mergeStrands=FALSE) {
  this <- extend(ProbeAffinityFile(...), "SnpProbeAffinityFile",
    mergeStrands=mergeStrands
  );

  # Parse attributes (all subclasses must call this in the constructor).
  setAttributesByTags(this)

  this;
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
setMethodS3("getCellIndices", "SnpProbeAffinityFile", function(this, ..., unlist=FALSE) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  cdfMergeStrands <- affxparser::cdfMergeStrands


  # Argument 'unlist':
  unlist <- Arguments$getLogical(unlist);

  # Supported case?
  mergeStrands <- this$mergeStrands;
  if (unlist && mergeStrands) {
    throw("Unsupported request: Argument 'unlist' have to be TRUE when parameter 'mergeStrands' is TRUE: ", unlist);
  }

  cells <- NextMethod("getCellIndices");

  # Merge strands?
  if (mergeStrands) {
    cells <- .applyCdfGroups(cells, cdfMergeStrands);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Unlist?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (unlist) {
    cells <- unlist(cells, use.names=FALSE);
  }


  cells;
}, protected=TRUE) # getCellIndices()


setMethodS3("setMergeStrands", "SnpProbeAffinityFile", function(this, status, ...) {
  this$mergeStrands <- status;
}, protected=TRUE)


############################################################################
# HISTORY:
# 2006-09-11
# o Created.
############################################################################
