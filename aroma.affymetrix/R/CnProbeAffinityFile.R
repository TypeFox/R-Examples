###########################################################################/**
# @RdocClass CnProbeAffinityFile
#
# @title "The CnProbeAffinityFile class"
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
#   \item{...}{Arguments passed to @see "SnpProbeAffinityFile".}
#   \item{combineAlleles}{If @FALSE, allele A and allele B are treated
#      seperately, otherwise together.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("CnProbeAffinityFile", function(..., combineAlleles=FALSE) {
  this <- extend(SnpProbeAffinityFile(...), "CnProbeAffinityFile",
    combineAlleles=combineAlleles
  );

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
#             of @see "SnpProbeAffinityFile".}
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
setMethodS3("getCellIndices", "CnProbeAffinityFile", function(this, ..., unlist=FALSE) {
  # Argument 'unlist':
  unlist <- Arguments$getLogical(unlist);


  # Supported case?
  combineAlleles <- this$combineAlleles;
  if (unlist && combineAlleles) {
    throw("Unsupported request: Argument 'unlist' have to be TRUE when parameter 'combineAlleles' is TRUE: ", unlist);
  }


  cells <- NextMethod("getCellIndices");

  # If combining alleles, return only every second group.
  # In order to improve readability we merge the names of alleles groups
  # combined, e.g. groups 'C' and 'G' become group 'CG'.
  if (combineAlleles) {
    cells <- .applyCdfGroups(cells, function(groups) {
      ngroups <- length(groups);
      odds <- seq(from=1L, to=ngroups, by=2L);
      names <- names(groups);
      groups <- groups[odds];
      if (ngroups >= 2L) {
        evens <- seq(from=2L, to=ngroups, by=2L);
        names <- paste(names[odds], names[evens], sep="");
        names(groups) <- names;
      }
      groups;
    })
  } # if (combineAlleles)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Unlist?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (unlist) {
    cells <- unlist(cells, use.names=FALSE);
  }

  cells;
}, protected=TRUE) # getCellIndices()


setMethodS3("setCombineAlleles", "CnProbeAffinityFile", function(this, status, ...) {
  this$combineAlleles <- status;
}, protected=TRUE)


############################################################################
# HISTORY:
# 2006-09-12
# o Updated.
# 2006-09-11
# o Created.
############################################################################
