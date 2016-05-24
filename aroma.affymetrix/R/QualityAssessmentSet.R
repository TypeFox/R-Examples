###########################################################################/**
# @RdocClass QualityAssessmentSet
#
# @title "The QualityAssessmentSet class"
#
# \description{
#  @classhierarchy
#
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to constructor of @see "AffymetrixCelSet".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "KS"
#*/###########################################################################
setConstructorS3("QualityAssessmentSet", function(...) {
  extend(AffymetrixCelSet(...), "QualityAssessmentSet");
})

##########################################################################
# HISTORY:
# 2007-12-10
# o Now a QualityAssessmentSet is a plain AffymetrixCelSet.
# o Removed argument 'tags' to constructor of QualityAssessmentSet,
#   because it was never used anywhere.
##########################################################################
