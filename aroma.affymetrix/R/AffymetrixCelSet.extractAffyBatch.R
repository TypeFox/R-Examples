###########################################################################/**
# @set "class=AffymetrixCelSet"
# @RdocMethod extractAffyBatch
# @alias extractAffyBatch.ChipEffectSet
# @alias extractAffyBatch
#
# @title "Extracts an in-memory AffyBatch object from the CEL set"
#
# \description{
#  @get "title".
#  Note that any modifications done to the extract object will \emph{not}
#  be reflected in the original CEL set.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Argument passed to \code{ReadAffy()}
#     (@see "affy::read.affybatch").}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns an @see "affy::AffyBatch-class" object.
# }
#
# \details{
#  Since the \pkg{affy} package is making use of special CDF environment
#  packages, this method will warn if the needed package is missing and
#  explain that \pkg{affy} will later try to download and install it
#  automatically.
# }
#
# @author "HB"
#
# \seealso{
#   Internally @see "affy::read.affybatch" is used to read the data.
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("extractAffyBatch", "AffymetrixCelSet", function(this, ..., verbose=FALSE) {
  requireNamespace("affy") || throw("Package not loaded: affy")
  cleancdfname <- affy::cleancdfname
  ReadAffy <- affy::ReadAffy


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  cdf <- getCdf(this);
  chipType <- getChipType(cdf, fullname=FALSE);
  cdfPkgName <- cleancdfname(chipType);
  suppressWarnings({
    .require <- require # To please R CMD check
    res <- .require(cdfPkgName, character.only=TRUE);
  });
  if (!res) {
    warning("CDF enviroment package '", cdfPkgName, "' not installed. The 'affy' package will later try to download it from Bioconductor and install it.");
  }

  filenames <- getPathnames(this);
  verbose && enter(verbose, "Creating AffyBatch from ", length(filenames), " CEL files");
  verbose && cat(verbose, "Filenames: ", paste(filenames, collapse=", "));
  sampleNames <- getFullNames(this);
  verbose && cat(verbose, "Sample names: ", paste(sampleNames, collapse=", "));

  # Sanity check
  dups <- sort(sampleNames[duplicated(sampleNames)]);
  if (length(dups) > 0) {
    throw(sprintf("Cannot load %s as an AffyBatch. Detected %d files that share the same sample names: %s", class(this)[1], length(dups)+length(unique(dups)), paste(unique(dups), collapse=", ")));
  }

  # Specify ReadAffy() of 'affy' to avoid conflicts with the one
  # in 'oligo'.
  read.affybatch <- affy::read.affybatch;
  ReadAffy <- affy::ReadAffy;
  res <- ReadAffy(filenames=filenames, sampleNames=sampleNames, ..., verbose=as.logical(verbose));

  verbose && exit(verbose);
  res;
}) # extractAffyBatch()


setMethodS3("extractAffyBatch", "ChipEffectSet", function(this, ...) {
  throw("Cannot extract AffyBatch from an ", class(this)[1], " object because it contains estimates that are summarized over sets of probes, whereas an AffyBatch should contain probe-level signals: ", getPath(this));
}, protected=TRUE)



############################################################################
# HISTORY:
# 2010-11-17
# o ROBUSTNESS: Now extractAffyBatch() for AffymetrixCelSet asserts that
#   the sample names are unique, which affy::ReadAffy() requires.
#   Moreover, the sample names are now the fullnames not just the names.
# 2010-09-06
# o ROBUSTNESS: Added extractAffyBatch() for ChipEffectSet that gives an
#   informative error message explaining why it doesn't make sense to do so.
# 2006-10-02
# o Created. A first small step toward an interface to Bioconductor.
############################################################################
