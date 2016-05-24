###########################################################################/**
# @set "class=AffymetrixCdfFile"
# @RdocMethod isSnpChip
#
# @title "Static method to check if a chip is a mapping (SNP) chip"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns @TRUE if the chip type refers to a SNP array, otherwise @FALSE.
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("isSnpChip", "AffymetrixCdfFile", function(this, ...) {
  chipType <- getChipType(this);

  # First some hardwired return values
  if (regexpr("^Mapping(10K|50K|250K)_.*$", chipType) != -1)
    return(TRUE);

  if (regexpr("^Cent(Hind|Xba).*$", chipType) != -1)
    return(TRUE);

  if (regexpr("^GenomeWideSNP_.*$", chipType) != -1)
    return(TRUE);

  if (regexpr("^Cyto.*Array$", chipType) != -1)
    return(TRUE);

  # Then, check for genotype units
  types <- getUnitTypes(this, ...);
  hasSnpUnits <- any(types == 2);

  hasSnpUnits;
}, private=TRUE)



###########################################################################/**
# @RdocMethod getSnpNames
#
# @title "Gets the names of the SNP units"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character @vector.
# }
#
# @author "HB"
#
# \seealso{
#   @seemethod "getCnNames".
#   Internally, @seemethod "getUnitTypes".
#   is used.
#
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getSnpNames", "AffymetrixCdfFile", function(this, ...) {
  types <- getUnitTypes(this, ...);
  units <- (types == 2);
  getUnitNames(this, units=units, ...);
}, private=TRUE)




###########################################################################/**
# @RdocMethod getCnNames
#
# @title "Gets the names of the CN units"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character @vector.
# }
#
# @author "HB"
#
# \seealso{
#   @seemethod "getSnpNames".
#   Internally, @seemethod "getUnitTypes".
#   is used.
#
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getCnNames", "AffymetrixCdfFile", function(this, ...) {
  types <- getUnitTypes(this, ...);
  units <- (types == 5);
  getUnitNames(this, units=units, ...);
}, private=TRUE)




###########################################################################/**
# @RdocMethod nbrOfSnps
#
# @title "Gets the number of SNPs"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Additional arguments passed to @seemethod "getSnpNames".}
# }
#
# \value{
#   Returns an @integer.
# }
#
# @author "HB"
#
# \seealso{
#   Internally, @seemethod "getSnpNames" is used to identify SNPs.
#
#   @seeclass
# }
#*/###########################################################################
setMethodS3("nbrOfSnps", "AffymetrixCdfFile", function(this, ...) {
  length(getSnpNames(this, ...));
}, private=TRUE)



############################################################################
# HISTORY:
# 2013-04-23
# o SPEEDUP: Added ^Cyto.*Array$ to isSnpChip() for AffymetrixCdfFile.
# 2008-05-10
# o ROBUSTNESS: Added backward compatibility for cases when the cached
#   results has sets$nonSNPs as a list.
# 2008-03-26
# o CLEAN UP: getAlleleProbePairs() of AffymetrixCdfFile would print *all*
#   identified non-SNP cells in the verbose output, instead of using str().
# 2008-02-27
# o Now getAlleleProbePairs() also returns element 'nonSNPs' (unless NULL).
# 2008-02-21
# o Now getAlleleProbePairs() only consider SNPs with 2 or 4 groups, because
#   at least one custom SNP chip we've seen a few SNPs with also 6 groups
#   (which turned out to all have the same direction).
# o GENERALIZED: Now getSnpNames(), getCnNames(), getAlleleProbePairs(),
#   getAlleleProbePairs2(), and isSnpChip() all infer unit type (SNP or CN)
#   from the CDF unit type and no longer from the unit names.
# 2007-09-14
# o Added getCnNames().
# o Updated isSnpChip() to recognize 5.0 and 6.0 chips.
# o Update regular expression for getSnpNames().
# 2007-08-16
# o Now getAlleleProbePairs() of AffymetrixCdfFile processes the CDF in
#   chunks in order to save memory.  Before the GenomeWideSNP_6 CDF would
#   consume 1.5-2.0GB RAM, but now it is using less than 500MB.
# 2007-06-11
# o BUG FIX: getAlleleProbePairs2() used non-existing object 'name' instead
#   of 'basepair'.  getAlleleProbePairs2() is currently not used anyway.
# 2006-09-15
# o Adopted to the new aroma.affymetrix structure.
# 2006-07-21
# o Added getAllelePairProbes().
# o Added getAllelePairProbesets().
# 2006-06-04
# o Added getAllelePairs().
# 2006-05-31
# o Added getSnpNames() and nbrOfSnps().
# 2006-05-30
# o Added static fromFile() which tries to call ditto of all subclasses.
# o Added static isSnpChip().
# 2006-03-30
# o Updated according to affxparser.
# 2006-03-27
# o Added detailed Rdoc comments to getRelativeAlleleSignals().
# 2006-03-24
# o Added references to DM articles and Affymetrix manuals.
# o Further speed up by improve rearrangement of CDF structure. Now a Hind
#   chip takes about 11-13 minutes instead.  11 minutes compared with
#   35 hours is 190 times faster.
# o After several speed improvements (also in affxparser), estimation of DM
#   rank scores now takes about 15-18 minutes for the 100K Hind chip.
#   The first draft took 30-35 hours(!) and yesterday 60-80 minutes.  Note,
#   the first draft was not "stupid" code; there is always room for
#   improvement.
# o Defined a local colSums() in getDmRankScores() specialized for matrices.
#   The overhead of the default colSums() is about 50%.
# 2006-03-23
# o Moved all SNP related methods into the new class AffymetrixSnpCelFile.
# o Added getRelativeAlleleSignals().  Note, it was designed to be used
#   with the 10K SNP chips.  These are designed so that there are equal
#   number of forward and reverse quartets with matching offsets in both
#   strands.  This is not the case for the 100K chips and above.
# o Created.
############################################################################
