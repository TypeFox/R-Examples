setConstructorS3("Transform", function(..., .reqSetClass="AffymetrixCelSet") {
  extend(AromaTransform(..., .reqSetClass=.reqSetClass), "Transform");
}, abstract=TRUE)


setMethodS3("getOutputFiles", "Transform", function(this, pattern=NULL, ...) {
  # Argument 'pattern':
  if (is.null(pattern)) {
    # Default filename pattern find non-private (no dot prefix) CEL files.
    pattern <- "^[^.].*[.](cel|CEL)$";
  } else {
    pattern <- Arguments$getRegularExpression(pattern=pattern);
  }

  NextMethod("getOutputFiles", pattern=pattern);
}, protected=TRUE)



###########################################################################/**
# @set "class=Transform"
# @RdocMethod getOutputDataSet
#
# @title "Gets the transformed data set"
#
# \description{
#  @get "title", if processed.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{force}{If @TRUE, any in-memory cached results are ignored.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns an @see "aroma.core::AromaMicroarrayDataSet" or @NULL.
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getOutputDataSet", "Transform", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Getting output data set for ", class(this)[1]);

  # Inherit the CDF from the input data set.
  ds <- getInputDataSet(this);
  cdf <- getCdf(ds);
  args <- list(generic="getOutputDataSet", this, ...,
               cdf=cdf, checkChipType=FALSE);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Inherit certain arguments from the input data set
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # AD HOC (not using OO), but setting these arguments does speed
  # up things. /HB 2007-09-17
  # Note, this is also done in Transform for now, such it is not really
  # needed here.  However, in case it will be removed from there it still
  # makes sense to have it here.
  if (inherits(ds, "CnChipEffectSet"))
    args$combineAlleles <- ds$combineAlleles;
  if (inherits(ds, "SnpChipEffectSet"))
    args$mergeStrands <- ds$mergeStrands;

  verbose && cat(verbose, "Calling NextMethod:");
  verbose && str(verbose, args);
  args$verbose <- less(verbose,1);

  res <- do.call(NextMethod, args);

  # Let the set update itself
  if (!is.null(res)) {
    update2(res, ..., verbose=less(verbose,1));
  }

  verbose && exit(verbose);

  res;
})


############################################################################
# HISTORY:
# 2012-10-14
# o CLEANUP: Removed obsolete getOutputDataSetOLD20090509() for Transform.
# 2009-05-23
# o Now getOutputDataSet() of Transform may return NULL if the output
#   data set is empty. Before it gave an error say update2() is not
#   applicable.
# 2009-05-09
# o Updated getOutputDataSet() of Transform to work with the updated
#   superclass AromaTransform.
# 2008-05-31
# o BUG FIX: The recent updates to getOutputFiles() did also find private
#   files.
# 2008-05-23
# o Transform now inherits from platform-independent AromaTransform.
# o Removed some dependencies to CDFs.
# 2007-12-08
# o getOutputDataSet() of Transform was updated to utilize the new 'cdf'
#   argument in static fromFiles() of AffymetrixCelSet.  This way the
#   default is not queried (in case it does not exist).
# 2007-09-18
# o Now getOutputDataSet() of Transform carry down certain arguments from
#   the input data set. This will speed up things.
# 2007-09-12
# o Now getOutputDataSet() of Transform passes down '...'static
#   fromFiles() of the AffymetrixCelSet class being setup.
# o Now isDone() of Transform throws an error if too many output files are
#   found.  Before it used to return FALSE.
# 2007-09-05
# o Added test against generating an output path that is the same as the
#   path of the input data set.
# 2007-06-25
# o BUG FIX: When getOutputDataSet() retrieved the output data set, the chip
#   type of the CEL files would be validated against the path name, also when
#   then CDF of the input set was overriden.  Now the output data set is
#   setup using 'checkChipType=FALSE'.  Thanks Mark Robinson for
#   troubleshooting this.
# 2007-03-24
# o BUG FIX: getPath() created the root path before trying to expand
#   Windows shortcuts.
# 2007-02-28
# o Now getOutputData() of Transform make sure to pass down the CDF too.
# 2007-01-14
# o Added a test for "unknown" (=unused) arguments to constructor.
# 2007-01-07
# o BUG FIX: getOutputFiles() would return "private" (prefix '.') files too.
#   This caused for instance FragmentLengthNormalization to return FALSE
#   for isDone() after being average, because one too many files was found.
# 2007-01-06
# o Renamed to Transform (from Preprocessing).
# 2006-12-20
# o Now isDone() returns FALSE if not all output files are there.  Before
#   an exception was thrown.  This modification allows you to for instance
#   remove a quantile normalized output file, and when reprocessing the
#   data set, only that file will be processed. I made this change after
#   one file was corrupted in a large data set and I did not want to have
#   to reprocess the whole data set.
# 2006-12-08
# o Renamed from PreProcessor.
# 2006-12-07
# o Created from QuantileNormalizer.R.
############################################################################
