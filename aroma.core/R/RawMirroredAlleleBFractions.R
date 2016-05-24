###########################################################################/**
# @RdocClass RawMirroredAlleleBFractions
#
# @title "The RawMirroredAlleleBFractions class"
#
# \description{
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "RawAlleleBFractions".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#*/###########################################################################
setConstructorS3("RawMirroredAlleleBFractions", function(...) {
  extend(RawAlleleBFractions(...), "RawMirroredAlleleBFractions");
})

setMethodS3("plot", "RawMirroredAlleleBFractions", function(x, ..., ylim=c(0,1/2)+c(-0.2,0.2), ylab="Mirrored Allele B fraction") {
  NextMethod("plot", ylim=ylim, ylab=ylab);
})

setMethodS3("extractRawMirroredAlleleBFractions", "default", abstract=TRUE);


setMethodS3("extractRawMirroredAlleleBFractions", "RawAlleleBFractions", function(this, ...) {
  beta <- getSignals(this);
  dh <- abs(beta - 1/2);
  res <- clone(this);
  res <- setSignals(res, dh);
  class(res) <- c("RawMirroredAlleleBFractions", class(res));
  res;
})




############################################################################
# HISTORY:
# 2009-05-17
# o Created.
############################################################################
