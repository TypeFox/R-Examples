###########################################################################/**
# @set "class=AffymetrixCelSet"
# @RdocMethod extractFeatureSet
# @alias extractFeatureSet
# @aliasmethod extractSnpFeatureSet
# @alias extractSnpFeatureSet
#
# @title "Extracts CEL signals an in-memory FeatureSet object"
#
# \description{
#  @get "title" from a @see "AffymetrixCelSet" object.
#  Note that any modifications done to the extracted object will \emph{not}
#  be reflected in the original CEL set.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Additional argument passed to @see "oligo::read.celfiles".}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns an @see "oligoClasses::FeatureSet-class" object.
# }
#
# @author "HB"
#
# \seealso{
#   Internally @see "oligo::read.celfiles" is used to read the data.
#   To read \emph{summarized} data as a @see "Biobase::ExpressionSet-class"
#   object, see \code{\link[aroma.affymetrix:extractExpressionSet.ChipEffectSet]{*extractExpressionSet}()}.
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("extractFeatureSet", "AffymetrixCelSet", function(this, ..., verbose=FALSE) {
  requireNamespace("oligo") || throw("Package not loaded: oligo");
  read.celfiles <- oligo::read.celfiles


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Reading ", class(this)[1L], "as a FeatureSet");
  verbose && cat(verbose, "Number of data files: ", length(this));
  verbose && cat(verbose, "Chip type: ", getChipType(this));

  pathnames <- getPathnames(this);
  verbose2 <- as.logical(verbose);
  res <- read.celfiles(pathnames, ..., verbose=verbose2);

  verbose && exit(verbose);

  res;
}) # extractFeatureSet()


setMethodS3("extractSnpFeatureSet", "AffymetrixCelSet", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Reading complete SnpFeatureSet");
  res <- extractFeatureSet(this, ...);

  # Sanity check
  if (!inherits(res, "SnpFeatureSet")) {
    throw("The read data is not of class 'SnpFeatureSet': ", class(res)[1L]);
  }

  verbose && exit(verbose);

  res;
}, protected=TRUE) # extractSnpFeatureSet()


############################################################################
# HISTORY:
# 2013-04-27
# o Documented extractFeatureSet() for AffymetrixCelSet.
# o Renamed extractSnpFeatureSet() to extractFeatureSet().  Keeping
#   old one for backward compatibility.
# 2009-10-16
# o Created.
############################################################################
