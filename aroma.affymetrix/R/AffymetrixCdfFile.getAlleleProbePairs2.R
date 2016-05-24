setMethodS3("getAlleleProbePairs2", "AffymetrixCdfFile", function(this, ..., verbose=FALSE) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  cdfGetGroups <- affxparser::cdfGetGroups

  # Look up base::apply(); '::' is expensive
  base_apply <- base::apply;


  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  verbose && enter(verbose, "Identifying the probes stratified by allele basepairs");
  cdfFile <- getPathname(this);

  # Identify all possible allele pairs
  verbose && enter(verbose, "Loading all possible allele basepairs");
  # Use only units that are SNPs
  types <- getUnitTypes(this, verbose=verbose);
  units <- which(types == 2);

  # Read group names for the SNPs
  groupNames <- .readCdfGroupNames(cdfFile, units=units);
  uGroupNames <- unique(groupNames);
  verbose && exit(verbose);

  uGroupNames0 <- lapply(uGroupNames, FUN=function(x) {
    x <- matrix(x, nrow=2)
    if (ncol(x) == 2) {
      # Take the complement bases for the reverse strand
      x[,2] <- c(A="T", C="G", G="C", T="A")[x[,2]];
    }
    x;

  })

  uBasepairs0 <- lapply(uGroupNames0, FUN=function(x) {
    base_apply(x, MARGIN=2, FUN=sort);
  })

  uBasepairs1 <- lapply(uBasepairs0, FUN=function(x) {
    base_apply(x, MARGIN=2, FUN=paste, collapse="");
  })

  # Get all unique allele basepairs
  uBasepairs <- sort(unique(unlist(uBasepairs1)));
  uBasepairs <- strsplit(uBasepairs, split="");

  # Create basepairs to group names map
  map <- vector("list", length(uBasepairs));
  names(map) <- uBasepairs;
  bpIdx <- vector("list", length(uBasepairs0));
  for (bp in uBasepairs) {
    for (kk in 1:length(bpIdx)) {
      set <- uBasepairs0[[kk]];
      bpIdx[[kk]] <- kk + which(bp == set)/10;
    }
    map[[bp]] <- unlist(bpIdx);
  }

  # Read all of the CDF file
  verbose && enter(verbose, "Loading cell indices for all probepairs");
  cdfAll <- .readCdfCellIndices(cdfFile, units=units, stratifyBy="pm");
  # Not needed anymore
  units <- NULL;
  verbose && exit(verbose);

  verbose && enter(verbose, "Stratifying by unique allele basepairs");
  probes <- vector("list", length(map));
  for (kk in 1:length(map)) {
    basepair <- names(map)[kk];
    verbose && enter(verbose, "Allele basepair ", basepair);

    bpIdx <- map[[kk]];
    gIdx <- as.integer(bpIdx);
    sIdx <- round(10*(bpIdx - gIdx));

    verbose && cat(verbose, "Located in ", length(unique(gIdx)), " group(s).");

    idx <- lapply(groupNames, FUN=identical, basepair);
    idx <- which(unlist(idx, use.names=FALSE));
    cdf <- cdfAll[idx];

    cdf0 <- vector("list", length=4);
    for (gg in 1:4) {
      cells <- .applyCdfGroups(cdf, cdfGetGroups, gg);
      cells <- unlist(cells, use.names=FALSE);
      cdf0[[gg]] <- cells;
    }
    probes[[kk]] <- cdf0;
    names(probes)[kk] <- basepair;

    verbose && exit(verbose);
  }
  verbose && exit(verbose);

  probes;
}, private=TRUE) # getAlleleProbePairs2()



############################################################################
# HISTORY:
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
