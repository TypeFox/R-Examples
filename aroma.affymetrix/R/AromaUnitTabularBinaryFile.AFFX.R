# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# BEGIN: Platform specific
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("importFrom", "AromaUnitTabularBinaryFile", function(this, src, ...) {
  if (inherits(src, "AffymetrixNetAffxCsvFile")) {
    importFromAffymetrixNetAffxCsvFile(this, src, ...);
  } else if (inherits(src, "DChipGenomeInformation")) {
    importFromDChipGenomeInformation(this, src, ...);
  } else if (inherits(src, "GenomeInformation")) {
    importFromGenomeInformation(this, src, ...);
  } else if (inherits(src, "AffymetrixTabularFile")) {
    importFromAffymetrixTabularFile(this, src, ...);
  } else if (inherits(src, "GenericTabularFile")) {
    importFromGenericTabularFile(this, src, ...);
  } else {
    throw("Do not know how to import from an src of class ", class(src)[1]);
  }
})



setMethodS3("getCdf", "AromaUnitTabularBinaryFile", function(this, ..., force=FALSE, .old=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  cdf <- this$.cdf;
  if (force || is.null(cdf)) {
    if (.old) {
      # Generate all possible fullname 'chipTypes' and search for the existance
      # of a CDF with the longest name *and* that have the same number of units.
      # This is AD HOC! We should really store the full chiptype of the
      # corresponding CDF in the internal file header. /HB 2007-12-10

      verbose && enter(verbose, "Searching for a match CDF");

      verbose && cat(verbose, "Filename: ", getFilename(this));
      name <- getName(this, ...);
      tags <- getTags(this, collapse=NULL, ...);

      validator <- function(cdf, ...) {
        (nbrOfUnits(cdf) == nbrOfUnits(this));
      }
      pathname <- findByCdf2(chipType=name, tags=tags, validator=validator,
                                                    verbose=less(verbose, 1));
      if (is.null(pathname)) {
        throw("Failed to locate a CDF for ", class(this)[1],
              " that have ", nbrOfUnits, " units: ", getFullName(this));
      }

      cdf <- AffymetrixCdfFile$fromFile(pathname);

      verbose && exit(verbose);
    } else {
      chipType <- getChipType(this);
      nbrOfUnits <- nbrOfUnits(this);
      cdf <- AffymetrixCdfFile$byChipType(chipType, nbrOfUnits=nbrOfUnits);
    }

    this$.cdf <- cdf;
  }

  cdf;
})



###########################################################################/**
# @set "class=AromaUnitTabularBinaryFile"
# @RdocMethod allocateFromCdf
#
# @title "Creates an AromaUnitTabularBinaryFile mapping to a given CDF"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{cdf}{The @see "AffymetrixCdfFile" used as a template and from
#      which the (full) chip type is taken.}
#   \item{path}{The path where to store the new file.}
#   \item{tags}{A @character @vector of optional tags appended to
#      the filename.}
#   \item{footer}{A nested named @list structure of additional attributes
#      that are saved in the file footer after the mandatory ones.}
#   \item{...}{Additional arguments passed to \code{allocate()} of
#      @see "aroma.core::AromaTabularBinaryFile".}
# }
#
# \value{
#  Returns a @see "aroma.core::AromaUnitTabularBinaryFile" object.
# }
#
# @author "HB"
#
# \seealso{
#   To update to file footer afterwards, see \code{writeFooter()}.
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("allocateFromCdf", "AromaUnitTabularBinaryFile", function(static, cdf, ...) {
  # Argument 'cdf':
  cdf <- Arguments$getInstanceOf(cdf, "AffymetrixCdfFile");

  allocateFromUnitNamesFile(static, unf=cdf, ...);
}, static=TRUE)




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# END: Platform specific
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


setMethodS3("importFromGenericTabularFile", "AromaUnitTabularBinaryFile", abstract=TRUE);


setMethodS3("importFromAffymetrixTabularFile", "AromaUnitTabularBinaryFile", function(this, src, ...) {
  # Argument 'src':
  src <- Arguments$getInstanceOf(src, "AffymetrixTabularFile");

  importFromGenomeInformation(this, src, ...);
});


setMethodS3("importFromAffymetrixNetAffxCsvFile", "AromaUnitTabularBinaryFile", abstract=TRUE, protected=TRUE);

setMethodS3("importFromDChipGenomeInformation", "AromaUnitTabularBinaryFile", function(this, src, ...) {
  # Argument 'src':
  src <- Arguments$getInstanceOf(src, "DChipGenomeInformation");

  importFromGenomeInformation(this, src, ...);
})


setMethodS3("importFromGenomeInformation", "AromaUnitTabularBinaryFile", abstract=TRUE, protected=TRUE);


############################################################################
# HISTORY:
# 2009-02-10
# o Added optional validation of number of units to byChipType().
# 2008-07-21
# o BUG FIX: byChipType() of AromaUnitTabularBinaryFile failed to locate
#   a valid tabular file if more than one was found and it was not the
#   last one that was matching.
# 2008-05-19
# o Added platform-independent allocateFromUnitNamesFile() which now also
#   writes footer attribute 'platform'.
# 2008-02-13
# o Added and updated Rdoc comments.
# 2008-01-19
# o Now AromaUnitTabularBinaryFile gets the chip type from the file footer.
# o ROBUSTNESS: Now fromChipType() of AromaUnitTabularBinaryFile validates
#   that the number of units in the located file match the number of units
#   in the CDF located using the same search parameters.
# 2007-12-10
# o Currently a AromaUnitTabularBinaryFile (e.g. AromaUgpFile) does not
#   contain information about the "fullname" chip type, but only the basic
#   chip-type name, e.g. we cannot infer the full chip-type name from
#   'GenomeWideSNP_5,Full,r2.ugp', but only 'GenomeWideSNP_5'. The fullname
#   should be the same as the full chip-type name of the CDF used to define
#   the the unit map, e.g. 'GenomeWideSNP_5,Full.CDF'.
#   We should add a header (or footer) field in the file format that
#   indicates the full chip type.
#   However, until that is done, the best we can do is to turn to the ad
#   hoc solution of scanning for the CDF with the longest matching fullname,
#   if both 'GenomeWideSNP_5,Full.CDF' and 'GenomeWideSNP_5.CDF' exists,
#   the we match the former to 'GenomeWideSNP_5,Full,r2.ugp'.  The fullname
#   chip type of the UGP is then full chip-type name of the CDF.  NOTE,
#   there is major drawback with this.  If the user deletes the "full" CDF,
#   the above approach would all of a sudden return a different full name!
# o Added clearCache().
# 2007-09-14
# o Renames createFromCdf() to allocateFromCdf().
# 2007-09-13
# o Created from AromaUflFile.R.
############################################################################
