###########################################################################/**
# @RdocClass ProbeAffinityFile
#
# @title "The ProbeAffinityFile class"
#
# \description{
#  @classhierarchy
#
#  This class represents estimates of probe affinities in probe-level models.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "ParameterCelFile".}
#   \item{probeModel}{The specific type of probe model.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#
# \seealso{
#   An object of this class is typically obtained through the
#   \code{getProbeAffinityFile()} method for the @see "ProbeLevelModel" class.
# }
#
#*/###########################################################################
setConstructorS3("ProbeAffinityFile", function(..., probeModel=c("pm", "mm", "pm-mm", "min1(pm-mm)", "pm+mm")) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'probeModel':
  probeModel <- match.arg(probeModel);

  extend(ParameterCelFile(...), "ProbeAffinityFile",
    "cached:.firstCells" = NULL,
    probeModel = probeModel
  )
})

setMethodS3("as.character", "ProbeAffinityFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character");
  s <- c(s, sprintf("Parameters: %s", getParametersAsString(this)));
  s;
}, protected=TRUE)


setMethodS3("getParameters", "ProbeAffinityFile", function(this, ...) {
  params <- NextMethod("getParameters");
  params$.firstCells <- NULL;
  params$probeModel <- this$probeModel;
  params;
}, protected=TRUE)



###########################################################################/**
# @RdocMethod getCellIndices
#
# @title "Retrieves tree list of cell indices for a set of units"
#
# \description{
#   @get "title" from the associated CDF.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Additional arguments passed to \code{getCellIndices()}
#             of @see "AffymetrixCdfFile".}
#  \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns a @list structure, where each element corresponds to a unit.
#   If argument \code{unlist=TRUE} is passed, an @integer @vector is returned.
# }
#
# \seealso{
#   @seeclass
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("getCellIndices", "ProbeAffinityFile", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Getting cell indices");
  verbose && cat(verbose, "Probe model: ", this$probeModel);

  stratifyBy <- switch(this$probeModel, "pm"="pm", "mm"="mm", "pm-mm"="pm", "min1(pm-mm)"="pm", "pm+mm"="pm");
  verbose && cat(verbose, "stratifyBy: ", stratifyBy);

  cdf <- getCdf(this);

  res <- getCellIndices(cdf, ..., stratifyBy=stratifyBy,
                                             verbose=less(verbose, 5));
  verbose && exit(verbose);

  res;
}, protected=TRUE) # getCellIndices()



setMethodS3("readUnits", "ProbeAffinityFile", function(this, units=NULL, cdf=NULL, ...) {
  if (is.null(cdf)) {
    cdf <- getCellIndices(this, units=units);
  }

  # Note that the actually call to the decoding is done in readUnits()
  # of the superclass.
  NextMethod("readUnits", cdf=cdf, readStdvs=TRUE, readPixels=TRUE);
});


setMethodS3("updateUnits", "ProbeAffinityFile", function(this, units=NULL, cdf=NULL, data, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Updating units of ", class(this)[1]);
  verbose && cat(verbose, "Units:");
  verbose && str(verbose, units);

  if (is.null(cdf)) {
    cdf <- getCellIndices(this, units=units, verbose=less(verbose, 1));
  }

  # Note that the actually call to the encoding is done in updateUnits()
  # of the superclass.
  res <- NextMethod("updateUnits", cdf=cdf, data=data, verbose=less(verbose, 1));

  verbose && exit(verbose);

  res;
}, private=TRUE)






############################################################################
# HISTORY:
# 2010-02-20
# o Added verbose output to updateUnits() of ProbeAffinityFile.
# 2007-05-09
# o Removed writeSpatial().
# 2007-01-03
# o Renamed constructor argument 'model' to 'probeModel'.
# 2006-09-11
# o Update read- and updateUnits() to make use of getCellIndices().
# o Added getCellIndices().
# 2006-08-26
# o Added writeSpatial().
# 2006-08-25
# o Added findUnitsTodo().
# o Added getFirstCellIndices(). Since reading all cell indices can take
#   a while it is cached in memory, but also on file (in case we restart).
# o Created from LiWongProbeAffinityFile.  The RMA version is almost
#   identical so I made this a superclass of both.
############################################################################
