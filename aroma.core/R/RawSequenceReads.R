###########################################################################/**
# @RdocClass RawSequenceReads
#
# @title "The RawSequenceReads class"
#
# \description{
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#   \item{x}{An @integer @vector of length J specifying the read positions.}
#   \item{y}{An (optional) @integer @vector of length J specifying the number of reads at each position. Default is one read per position.}
#   \item{...}{Arguments passed to @see "RawGenomicSignals".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#*/########################################################################### 
setConstructorS3("RawSequenceReads", function(x=NULL, y=rep(1L, length(x)), ...) {
  extend(RawGenomicSignals(y=y, x=x, ...), "RawSequenceReads")
})

setMethodS3("nbrOfReads", "RawSequenceReads", function(this, ...) {
  nbrOfLoci(this, ...);
})

setMethodS3("binnedSums", "RawSequenceReads", function(this, ...) {
  binnedSmoothing(this, ..., FUN=sum);
})


setMethodS3("plot", "RawSequenceReads", function(x, ..., ylim=c(0,10), ylab="Reads") {
  NextMethod("plot", ylim=ylim, ylab=ylab);
})



############################################################################
# HISTORY:
# 2010-09-11
# o Added an explicit 'ylim' argument to plot() for RawSequenceReads. 
# 2009-07-02
# o Created from RawCopyNumbers.R.
############################################################################
