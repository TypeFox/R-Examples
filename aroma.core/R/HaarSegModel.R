###########################################################################/**
# @RdocClass HaarSegModel
#
# @title "The HaarSegModel class"
#
# \description{
#  @classhierarchy
#
#  This class represents the Haar wavelet-based segmentation (HaarSeg)
#  model [1].
# }
#
# @synopsis
#
# \arguments{
#   \item{cesTuple}{A @see "CopyNumberDataSetTuple".}
#   \item{breaksFdrQ}{Default tuning parameters specific to the HaarSeg
#         algorithm.}
#   \item{...}{Arguments passed to the constructor of
#              @see "CopyNumberSegmentationModel".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#
# \seealso{
#  @see "CopyNumberSegmentationModel".
# }
#
# \references{
#   [1] Ben-Yaacov E. and Eldar YC. \emph{A fast and flexible method for the segmentation of aCGH data}, Bioinformatics, 2008.
#   \url{http://www.ee.technion.ac.il/Sites/People/YoninaEldar/Info/software/HaarSeg.htm}
# }
#*/###########################################################################
setConstructorS3("HaarSegModel", function(cesTuple=NULL, ..., breaksFdrQ=0.0001) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Load required packages
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(cesTuple)) {
    # Early error, iff package is not installed
    requireNamespace("HaarSeg") || throw("Package not loaded: HaarSeg");
  }

  # Argument 'breaksFdrQ':
  breaksFdrQ <- Arguments$getDouble(breaksFdrQ, range=c(0,1));

  extend(CopyNumberSegmentationModel(cesTuple=cesTuple, ...), "HaarSegModel",
    .breaksFdrQ = breaksFdrQ
  )
})


setMethodS3("getAsteriskTags", "HaarSegModel", function(this, ...) {
  NextMethod("getAsteriskTags", tag="HAAR");
}, protected=TRUE)


setMethodS3("getFitFunction", "HaarSegModel", function(this, ...) {
  segmentByHaarSeg;
}, protected=TRUE)


##############################################################################
# HISTORY:
# 2009-05-16
# o Added getFitFunction().  Removed fitOne().
# 2008-05-14
# o Moved drawCnRegions(), extractCopyNumberRegions() and
#   extractRawCopyNumbers() for HaarSeg to aroma.core v1.0.6 (will
#   eventually end up in aroma.cn).
# 2008-01-26
# o Reordered constructor arguments.
# 2008-12-31
# o Removing non-finite data points before passing to haarSeg().
# 2008-12-17
# o Now using the HaarSeg package (put together by HB).
# 2008-12-16
# o Created.
##############################################################################
