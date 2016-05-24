###########################################################################/**
# @RdocClass ProbeLevelTransform
#
# @title "The ProbeLevelTransform class"
#
# \description{
#  @classhierarchy
#
#  This abstract class represents a transformation methods that transforms
#  probe-level signals, typically intensities.
# }
#
# @synopsis
#
# \arguments{
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
setConstructorS3("ProbeLevelTransform", function(...) {
  extend(Transform(...), "ProbeLevelTransform")
}, abstract=TRUE)


setMethodS3("getRootPath", "ProbeLevelTransform", function(this, ...) {
  # Ad hoc fix: /HB 2007-04-11
  ds <- getInputDataSet(this);
  if (inherits(ds, "ChipEffectSet"))
    return("plmData");

  "probeData";
}, protected=TRUE)



############################################################################
# HISTORY:
# 2006-12-08
# o Created.
############################################################################
