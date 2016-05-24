###########################################################################/**
# @set "class=matrix"
# @RdocMethod fitCalMaTeCNprobes
# @alias fitCalMaTeCNprobes
#
# @title "Normalizes non-polymorphic copy number loci according to the CalMaTe method"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{dataT}{A JxI @numeric @matrix, where J is the number of loci
#                      and I is the number of samples.}
#  \item{references}{A @integer @vector with elements in [1,I] specifying
#     which samples should be used as the reference set.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @numeric @vector of length J.
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("fitCalMaTeCNprobes", "matrix", function(dataT, references, ...) {
  # This is an internal function. Because of this, we will assume that
  # all arguments are valid and correct.  No validation will be done.

  # Extract the reference samples
  Tref <- dataT[,references, drop=FALSE];
  
  TR <- rowMedians(Tref, na.rm=TRUE);

  res <- 2 * dataT/TR;
  
  res;
}, protected=TRUE) # fitCalMaTeCNprobes()


###########################################################################
# HISTORY:
# 2011-11-29 [MO]
# o Change matrix "T" by "dataT".
# o Clear code, commented lines removed.
# 2010-08-02 [HB]
# o ROBUSTNESS: fitCalMaTeCNprobes() can now also handle missing values.
# o Made into an S3 method for matrix:es.
# 2010-06-22 [MO]
# o Created.
###########################################################################
