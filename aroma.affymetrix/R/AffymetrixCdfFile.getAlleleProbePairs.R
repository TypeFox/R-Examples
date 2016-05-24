###########################################################################/**
# @set "class=AffymetrixCdfFile"
# @RdocMethod getAlleleProbePairs
#
# @title "Gets the indices of probepairs with the same pair of SNP nucleotides"
#
# \description{
#   @get "title".
#   Note that the order of allele A and allele B is irrelevant.
#   For instance, all probepairs with nucleotides (A,T) are calibrated
#   together with all probepairs with nucleotides (T,A) reversed.
# }
#
# @synopsis
#
# \arguments{
#  \item{verbose}{A @logical or a @see "R.utils::Verbose" object.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a named @list where each element is a two-column @matrix where
#   the column names are the nucleotides for the two alleles.
# }
#
# \section{Benchmarking}{
#   On an IBM Thinkpad A31 with 1.8GHz and 1GB RAM:
#   \itemize{
#    \item{Mapping10K\_Xba142}{10208 units & 432964 cells: 11 seconds.}
#    \item{Mapping50K\_Xba240}{58960 SNPs & 589,600 (PMA,PMB) probe pairs: 11 seconds.}
#   }
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getAlleleProbePairs", "AffymetrixCdfFile", function(this, units=NULL, ignoreOrder=TRUE, force=FALSE, verbose=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Identifying the probes stratified by allele basepairs");
  on.exit(verbose && exit(verbose));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached results?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  chipType <- getChipType(this);
  key <- list(method="getAlleleProbePairs", class=class(this)[1], version="2008-02-27", chipType=chipType, units=units, ignoreOrder=ignoreOrder);
  if (getOption(aromaSettings, "devel/useCacheKeyInterface", FALSE)) {
    key <- getCacheKey(this, method="getAlleleProbePairs", chipType=chipType, units=units, ignoreOrder=ignoreOrder);
  }
  dirs <- c("aroma.affymetrix", chipType);
  if (!force) {
    probeSets <- loadCache(key=key, dirs=dirs);
    if (!is.null(probeSets)) {
      # Backward compatibility; remove October 2008. /HB
      if (is.list(probeSets$nonSNPs))
        break;
      verbose && cat(verbose, "Loaded from file cache");
      gc <- gc();
      verbose && print(verbose, gc);
      return(probeSets);
    }
  }

  cdfFile <- getPathname(this);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify all possible allele pairs
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Loading all possible allele basepairs");

  verbose && enter(verbose, "Identifying compatible SNP units");
  # Use only units that are SNPs...
  types <- getUnitTypes(this, verbose=less(verbose, 1));
  unitsAll <- which(types == 2);
  verbose && cat(verbose, "Number of SNP units: ", length(unitsAll));

  # ...and with either 2 or 4 groups
  unitSizes <- nbrOfGroupsPerUnit(this, units=unitsAll);
  verbose && cat(verbose, "Detected unit sizes:");
  verbose && print(verbose, table(unitSizes));

  unitsAll <- unitsAll[unitSizes %in% c(2,4)];
  # Not needed anymore
  unitSizes <- NULL;
  verbose && cat(verbose, "Number of SNP units with 2 or 4 groups: ",
                                                        length(unitsAll));

  verbose && exit(verbose);

  gc <- gc();

  # Operate only on a subset of probes?
  if (!is.null(units)) {
    unitsAll <- intersect(unitsAll, units);
  }
  units <- unitsAll;
  # Not needed anymore
  unitsAll <- NULL;
  nunits <- length(units);

  verbose && cat(verbose, "Number of SNP units to query: ", nunits);
  if (nunits == 0)
    return(NULL);

  # Read group names for these SNPs
  verbose && enter(verbose, "Retrieving group names");
  groupNames <- .readCdfGroupNames(cdfFile, units=units);
  # Save memory by removing names. [55Mb -> 44Mb]
  names(groupNames) <- NULL;
  # Save memory by converting to integers. [44Mb -> 11Mb]
  levels <- as.integer(1:4);
  names(levels) <- c("A", "C", "G", "T");
  groupNames <- lapply(groupNames, FUN=function(s) {
    s <- levels[s];
    names(s) <- NULL;
    s;
  });

  uGroupNames <- unique(groupNames);
  # Order by basepairs so that the verbose output is easier to read
  o <- order(as.integer(sapply(uGroupNames, FUN=paste, collapse="")));
  uGroupNames <- uGroupNames[o];
  # Not needed anymore
  o <- NULL;

  gc <- gc();
  verbose && print(verbose, gc);

  verbose && cat(verbose, "Unique group names:");
  verbose && str(verbose, lapply(uGroupNames, FUN=function(x) names(levels[x])), vec.len=8);

  verbose && exit(verbose);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read all of the CDF file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Loading cell indices for probepairs for requested units");
  nbrOfUnitsPerChunk <- 100e3;
#  nbrOfUnitsPerChunk <- 6000;
  nunits <- length(units);
  nbrOfChunks <- ceiling(nunits / nbrOfUnitsPerChunk);
  uu <- 1:nbrOfUnitsPerChunk;
  unitsTodo <- units;
  count <- 1;
  cells0 <- list();
  cdfAll <- list();
  while (length(unitsTodo) > 0) {
    verbose && enter(verbose, sprintf("Chunk #%d of %d", count, nbrOfChunks));
    if (length(unitsTodo) <= nbrOfUnitsPerChunk)
      uu <- 1:length(unitsTodo);

    verbose && cat(verbose, "Units: ");
    verbose && str(verbose, unitsTodo[uu]);

    cdfAll0 <- .readCdfCellIndices(cdfFile, units=unitsTodo[uu], stratifyBy="pm");
    unitsTodo <- unitsTodo[-uu];

    # Save memory by removing names. [309Mb -> 298Mb]
    names(cdfAll0) <- NULL;

    cells0[[count]] <- unlist(cdfAll0, use.names=FALSE);

    # Save memory by flattening structure. [298Mb -> 51Mb(!)]
    # TODO: Add support to do this already in affxparser?! /HB 2006-07-22
    cdfAll0 <- lapply(cdfAll0, FUN=function(unit) {
      groups <- .subset2(unit, 1);
      names(groups) <- NULL;
      lapply(groups, FUN=.subset2, 1);
    });

    gc <- gc();

    cdfAll <- c(cdfAll, cdfAll0);

    # Not needed anymore
    cdfAll0 <- NULL;
    gc <- gc();
    verbose && print(verbose, gc);

    count <- count + 1;
    verbose && exit(verbose);
  } # while(...)

  # Not needed anymore
  # Not needed anymore
  units <- unitsTodo <- uu <- NULL;
  gc <- gc();
  verbose && print(verbose, gc);

  cells0 <- unlist(cells0, use.names=FALSE);
  gc <- gc();
  cells0 <- sort(cells0);
  gc <- gc();

  nbrOfCells <- length(cells0);
  verbose && printf(verbose, "Identified %d (PM_A,PM_B) pairs in %d units, i.e. on average %.2g probe pairs/units\n", round(nbrOfCells/2), nunits, (nbrOfCells/2)/nunits);

  if (length(cdfAll) != nunits) {
    throw("Internal error: Expected ", nunits, " units, but got ", length(cdfAll));
  }

  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Group all units with the same allele basepairs
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Stratifying by unique allele basepairs");
  probeSets <- vector("list", length(uGroupNames));
  for (kk in 1:length(uGroupNames)) {
    name <- uGroupNames[[kk]];
    basepair <- paste(names(levels)[name[1:2]], collapse="");
    verbose && enter(verbose, sprintf("Allele basepair %s (%d of %d)", basepair, kk, length(uGroupNames)));

    idx <- sapply(groupNames, FUN=identical, name);
    idx <- which(idx);
    if (verbose) {
      bpNames <- matrix(names(levels)[name], nrow=2);
      bpNames <- paste(bpNames[1,], bpNames[2,], sep="");
      verbose && cat(verbose, "Allele pairs: ", paste(bpNames, collapse=","));
      # Not needed anymore
      bpNames <- NULL;
      verbose && cat(verbose, "Number of units: ", length(idx));
    }
    cdf <- cdfAll[idx];
    cdfAll[idx] <- NA;  # Not needed anymore (save memory)
    # Not needed anymore
    idx <- NULL;
#    gc <- gc();

    cdf0 <- vector("list", length=length(name));
    for (gg in 1:length(name)) {
      cdf0[[gg]] <- unlist(lapply(cdf, FUN=.subset2, gg), use.names=FALSE);
    }

    # Not needed anymore
    cdf <- NULL;
    probeSets[[kk]] <- cdf0;
    # Not needed anymore
    cdf0 <- NULL;
    names(probeSets)[kk] <- basepair;

#    gc <- gc();
#    verbose && print(verbose, gc);

    verbose && exit(verbose);
  }
  # Not needed anymore
  cdfAll <- NULL;
  gc <- gc();
  verbose && print(verbose, gc);
  verbose && exit(verbose);

  # Assert correctness
  verbose && enter(verbose, "Asserting correctness part I", level=-20);
  nbrOfCells2 <- length(unlist(probeSets, use.names=FALSE));
  if (nbrOfCells2 != nbrOfCells) {
    throw("Internal error: Excepted ", nbrOfCells, " indices: ", nbrOfCells2);
  }
  if (!identical(sort(unlist(probeSets, use.names=FALSE)), cells0)) {
    throw("Internal error: Mismatching probes.");
  }
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify equivalent groups
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Putting equivalent groups together");

  probeSets2 <- list();
  for (kk in 1:length(probeSets)) {
    bp <- names(probeSets)[kk];
    value <- probeSets[[kk]];

    nbrOfPairs <- length(value)/2;

    # The below assumes that every 2nd group pair is "reversed".
    # Should really make use of the group 'direction' in the CDF,
    # but that is really slow.  For the CDFs we've checked, all
    # 2 and 4 groups follow this assumption. /HB 2008-02-21
    for (ll in seq_len(nbrOfPairs)) {
      value2 <- probeSets2[[bp]];
      if (is.null(value2))
        value2 <- vector("list", length=2);
      value2[[1]] <- c(value2[[1]], value[[1]]);
      value2[[2]] <- c(value2[[2]], value[[2]]);
      probeSets2[[bp]] <- value2;
      # Not needed anymore
      value2 <- NULL;
      bp <- strsplit(bp, split="")[[1]];
      bp <- c(A="T", C="G", G="C", T="A")[bp];
      bp <- paste(bp, collapse="");
      value <- value[-(1:2)];
    }
  }
  verbose && cat(verbose, "Probe pairs: ",
                               paste(sort(names(probeSets2)), collapse=", "));
  verbose && exit(verbose);

  # Assert correctness
  verbose && enter(verbose, "Asserting correctness part II", level=-20);
  nbrOfCells2 <- length(unlist(probeSets, use.names=FALSE));
  if (nbrOfCells2 != nbrOfCells) {
    throw("Internal error: Excepted ", nbrOfCells, " indices: ", nbrOfCells2);
  }
  if (!identical(sort(unlist(probeSets, use.names=FALSE)), cells0)) {
    throw("Internal error: Mismatching probes.");
  }
  verbose && exit(verbose);

  if (ignoreOrder) {
    verbose && enter(verbose, "Putting AB and BA groups together");
    # Not needed anymore
    probeSets <- NULL;
    gc <- gc();
    pairs <- strsplit(names(probeSets2), split="");
    pairs <- lapply(pairs, FUN=function(x) paste(sort(x), collapse=""));
    pairs <- unlist(pairs);
    uPairs <- sort(unique(pairs));
    verbose && cat(verbose, "Probe pairs (ignoring order): ",
                                                paste(uPairs, collapse=", "));
    probeSets <- list();
    for (pair in uPairs) {
      idx <- which(pairs == pair);
      basepairs <- sort(names(probeSets2)[idx]);
      probeSets[[pair]] <- probeSets2[basepairs];
    }
    # Not needed anymore
    probeSets2 <- NULL;
    verbose && exit(verbose);

    verbose && enter(verbose, "Combining AB and BA groups");
    # Join AB with BA.
    for (kk in 1:length(probeSets)) {
      values <- probeSets[[kk]];
      if (length(values) == 1) {
        values <- values[[1]];
      } else {
        values[[1]][[1]] <- c(values[[1]][[1]], values[[2]][[2]]);
        values[[1]][[2]] <- c(values[[1]][[2]], values[[2]][[1]]);
        values <- values[[1]];
      }
      probeSets[[kk]] <- values;
    }
    # Not needed anymore
    values <- NULL;
    verbose && exit(verbose);
  } else {
    probeSets <- probeSets2;
    # Not needed anymore
    probeSets2 <- NULL;
  }

  # Assert correctness
  verbose && enter(verbose, "Asserting correctness part III", level=-20);
  nbrOfCells2 <- length(unlist(probeSets, use.names=FALSE));
  if (nbrOfCells2 != nbrOfCells) {
    throw("Internal error: Excepted ", nbrOfCells, " indices: ", nbrOfCells2);
  }
  if (!identical(sort(unlist(probeSets, use.names=FALSE)), cells0)) {
    throw("Internal error: Mismatching probes.");
  }
  gc <- gc();
  verbose && exit(verbose);

  verbose && enter(verbose, "Reformatting to matrices");
  # Order indices by allele A (just for beauty)
  for (kk in 1:length(probeSets)) {
    verbose && enter(verbose, sprintf("Group #%d of %d", kk, length(probeSets)));
    values <- probeSets[[kk]];
    values <- matrix(c(values[[1]], values[[2]]), ncol=2);
    colnames(values) <- strsplit(names(probeSets)[kk], split="")[[1]];
    o <- order(values[,1]);
    values <- values[o,];
    probeSets[[kk]] <- values;
    verbose && exit(verbose);
  }
  # Not needed anymore
  values <- o <- NULL;
  gc <- gc();
  if (isVisible(verbose, level=-20))
    verbose && str(verbose, probeSets, level=-20);
  verbose && exit(verbose);

  # Assert correctness
  verbose && enter(verbose, "Asserting correctness part IV", level=-20);
  nbrOfCells2 <- length(unlist(probeSets, use.names=FALSE));
  if (nbrOfCells2 != nbrOfCells) {
    throw("Internal error4: Excepted ", nbrOfCells, " indices: ", nbrOfCells2);
  }
  if (!identical(sort(unlist(probeSets, use.names=FALSE)), cells0)) {
    throw("Internal error: The identified set of indices for various allele probe pairs does not match the original set of cell indices.");
  }
  verbose && exit(verbose);


  verbose && enter(verbose, "Identifying indices for all non-SNP PM cells");
## OLD WAY!
##        unitNames <- getUnitNames(this);
##        snpNames <- getSnpNames(this);
##        nonSnpUnits <- which(!(unitNames %in% snpNames));
##        # Not needed anymore
##        unitNames <- snpNames <- NULL;

  # Identifying all units types
  unitTypes <- getUnitTypes(this, verbose=less(verbose,1)); # Takes time
  verbose && cat(verbose, "Table of identified unit types:");
  verbose && print(verbose, table(unitTypes));
  nonSnpUnits <- which(unitTypes != 2);  # '2 == genotype unit'

  if (length(nonSnpUnits) > 0) {
    cells <- getCellIndices(this, units=nonSnpUnits,
                useNames=FALSE, unlist=TRUE, verbose=less(verbose,1));
  } else {
    cells <- NULL;
  }
  verbose && cat(verbose, "Identified non-SNP units:");
  verbose && str(verbose, cells);
  probeSets$nonSNPs <- cells;
  # Not needed anymore
  cells <- NULL;
  verbose && exit(verbose);


  # Save cache to file
  comment <- key[c("method", "class", "chipType")];
  comment <- paste(names(comment), comment, sep="=");
  comment <- paste(comment, collapse=", ");
  saveCache(probeSets, key=key, comment=comment, dirs=dirs);

  probeSets;
}, private=TRUE) # getAlleleProbePairs()




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
