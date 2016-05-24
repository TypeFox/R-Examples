###########################################################################/**
# @RdocClass SegmentedCopyNumbers
#
# @title "The SegmentedCopyNumbers class"
#
# \description{
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "RawCopyNumbers".}
#   \item{states}{A @function returning the copy-number states given a
#     @vector of locus positions.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @examples "../incl/SegmentedCopyNumbers.Rex"
#
# @author
#*/########################################################################### 
setConstructorS3("SegmentedCopyNumbers", function(..., states=NULL) {
  this <- extend(RawCopyNumbers(...), c("SegmentedCopyNumbers", 
                                   uses("SegmentedGenomicSignalsInterface")));
  this <- setStates(this, states=states);
  this;
})

############################################################################
# HISTORY:
# 2012-03-01
# o No longer assumes that class provides reference variables.
# 2009-06-10
# o Now SegmentedCopyNumbers implements the SegmentedGenomicSignalsInterface
#   which makes this class a light-weight class.
# 2009-05-16
# o Now all methods of SegmentedCopyNumbers() coerce numerics only if
#   necessary, i.e. it keeps integers if integers, otherwise to doubles.
#   This is a general design of aroma.* that saves some memory.
# 2009-04-06
# o Now binnedSmoothingByState() of SegmentedCopyNumbers uses 
#   extractSubsetByState() and then binnedSmoothing() on that object.  
#   This makes the code slightly less redundant.
# 2009-02-19
# o Adopted to make use of new RawGenomicSignals.
# 2009-02-16
# o Now getStates() also passes the optional 'name' field to the "truth"
#   function.
# 2009-02-08
# o Added getUniqueStates().
# 2009-02-07
# o Added extractSubsetByState().
# o Created.
############################################################################
