###########################################################################/**
# @RdocClass ChipEffectTransform
#
# @title "The ChipEffectTransform class"
#
# \description{
#  @classhierarchy
#
#  This abstract class represents a transform that transforms chip-effect
#  estimates obtained from probe-level modelling.
# }
#
# @synopsis
#
# \arguments{
#   \item{dataSet}{The input data set as an @see "ChipEffectSet".}
#   \item{...}{Arguments passed to the constructor of @see "Transform".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \details{
#   Subclasses must implement the \code{process()} method.
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("ChipEffectTransform", function(dataSet=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  if (!is.null(dataSet)) {
    dataSet <- Arguments$getInstanceOf(dataSet, "ChipEffectSet");
  }

  extend(Transform(dataSet=dataSet, ...), "ChipEffectTransform")
}, abstract=TRUE)


setMethodS3("getRootPath", "ChipEffectTransform", function(this, ...) {
  "plmData";
}, protected=TRUE)


############################################################################
# HISTORY:
# 2009-05-09
# o CLEANUP: Removed getOutputDataSet() of ChipEffectTransform.  It is now
#   taken care of by the superclass.
# 2007-09-18
# o Now getOutputDataSet() of Transform carry down certain arguments from
#   the input data set. This will speed up things.
# 2006-12-08
# o Created.
############################################################################
