###########################################################################/**
# @RdocClass RawAlleleBFractions
#
# @title "The RawAlleleBFractions class"
#
# \description{
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "RawGenomicSignals".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#*/########################################################################### 
setConstructorS3("RawAlleleBFractions", function(...) {
  extend(RawGenomicSignals(...), "RawAlleleBFractions")
})

setMethodS3("plot", "RawAlleleBFractions", function(x, ..., ylim=c(0,1)+c(-0.2,0.2), ylab="Allele B fraction") {
  NextMethod("plot", ylim=ylim, ylab=ylab);
})

setMethodS3("extractRawAlleleBFractions", "default", abstract=TRUE);



############################################################################
# HISTORY:
# 2009-05-10
# o Created from RawCopyNumbers.R.
############################################################################
