###########################################################################/**
# @RdocClass BackgroundCorrection
#
# @title "The BackgroundCorrection class"
#
# \description{
#  @classhierarchy
#
#  This class represents a background adjustment function.
#
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to the constructor of
#     @see "ProbeLevelTransform".}
#   \item{subsetToUpdate}{The probes to be updated.
#     If @NULL, all probes are updated.}
#   \item{typesToUpdate}{Types of probes to be updated.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "KS, HB"
#*/###########################################################################
setConstructorS3("BackgroundCorrection", function(..., subsetToUpdate=NULL, typesToUpdate=NULL) {
  extend(ProbeLevelTransform(...), "BackgroundCorrection",
    .subsetToUpdate = subsetToUpdate,
    .typesToUpdate = typesToUpdate
  )
})


setMethodS3("getSubsetToUpdate", "BackgroundCorrection", function(this, ...) {
  this$.subsetToUpdate;
}, private=TRUE)


setMethodS3("getParameters", "BackgroundCorrection", function(this, ...) {
  # Get parameters from super class
  params <- NextMethod("getParameters");

  # Get parameters of this class
  params2 <- list(
    subsetToUpdate = this$.subsetToUpdate,
    typesToUpdate = this$.typesToUpdate
  );

  # Append the two sets
  params <- c(params, params2);

  params;
}, protected=TRUE)


###########################################################################/**
# @RdocMethod process
#
# @title "Processes the data set"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{force}{If @TRUE, data already processed is re-processed,
#       otherwise not.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @double @vector.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("process", "BackgroundCorrection", abstract=TRUE);



############################################################################
# HISTORY:
# 2007-03-21
# o Created.
############################################################################
