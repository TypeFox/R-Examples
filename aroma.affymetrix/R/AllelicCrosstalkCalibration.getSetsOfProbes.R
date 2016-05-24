setMethodS3("getSetsOfProbes", "AllelicCrosstalkCalibration", function(this, ..., version=c(0,1,3,4), fakeSymmetry=FALSE, force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'version':
  version <- Arguments$getIntegers(version);
  version <- version[1];
  if (version == 0) {
    if (this$.pairBy == "CDF") {
      version <- 1;
    } else if (this$.pairBy == "sequence") {
      version <- 4;
    } else {
      version <- 4;
    }
  }

  # Argument 'fakeSymmetry':
  fakeSymmetry <- Arguments$getLogical(fakeSymmetry);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  setsOfProbes <- this$.setsOfProbes;

  if (force || is.null(setsOfProbes)) {
    verbose && enter(verbose, "Identifying sets of pairs of cell indices");
    dataSet <- getInputDataSet(this);
    cdf <- getCdf(dataSet);
    chipType <- getChipType(cdf);
    verbose && cat(verbose, "Chip type: ", chipType);

    params <- getParameters(this);
    mergeShifts <- params$mergeShifts;
    verbose && cat(verbose, "Merge shifts: ", mergeShifts);
    B <- params$B;
    verbose && cat(verbose, "Number of nucleotides: ", B);

    # Assert that needed annotation files exist
    needAcs <- FALSE;
    if (version == 4) {
      needAcs <- TRUE;
    } else if (B > 1 || !mergeShifts) {
      needAcs <- TRUE;
    }

    if (needAcs) {
      verbose && enter(verbose, "Locating AromaCellSequenceFile");
      chipType <- getChipType(cdf, fullname=FALSE);
      nbrOfCells <- nbrOfCells(cdf);
      verbose && cat(verbose, "Chip type: ", chipType);
      verbose && cat(verbose, "Number of cells: ", nbrOfCells);
      acs <- AromaCellSequenceFile$byChipType(chipType, nbrOfCells=nbrOfCells);
      verbose && print(verbose, acs);
      verbose && exit(verbose);
    } else {
      acs <- NULL;
    }


    if (version == 1) {
      setsOfProbes <- getAlleleProbePairs(cdf, verbose=verbose);
      # Coerce to new structure. /HB 2008-09-02
      nonSNPs <- setsOfProbes$nonSNPs;
      setsOfProbes$nonSNPs <- NULL;
      setsOfProbes <- lapply(setsOfProbes, FUN=t);
      pairs <- names(setsOfProbes);
      pairs <- strsplit(pairs, split="", fixed=TRUE);
      pairs <- sapply(pairs, FUN=paste, collapse="/");
      names(setsOfProbes) <- pairs;
      setsOfProbes <- list(snps=setsOfProbes, nonSNPs=nonSNPs);
      # Not needed anymore
      nonSNPs <- pairs <- NULL;
    } else if (version == 3) {
      setsOfProbes <- getAlleleProbePairs3(cdf, ignoreOrder=TRUE, verbose=verbose);
      snps <- setsOfProbes$snps;
      for (kk in seq_along(snps)) {
        cells <- snps[[kk]][3:4,,drop=FALSE];
        o <- order(cells[1,]); 
        cells <- cells[,o,drop=FALSE];
        snps[[kk]] <- cells;
        # Not needed anymore
        cells <- o <- NULL;
      }
      pairs <- names(snps);
      pairs <- strsplit(pairs, split="", fixed=TRUE);
      pairs <- sapply(pairs, FUN=paste, collapse="/");
      names(snps) <- pairs;
      setsOfProbes$snps <- snps;
      # Not needed anymore
      snps <- pairs <- NULL;
    } else if (version == 4) {
      verbose && enter(verbose, "Identifying cell indices for all non-SNP units");
      unitTypes <- getUnitTypes(cdf, verbose=verbose);
      units <- which(unitTypes != 2);
      # Not needed anymore
      unitTypes <- NULL; 
      verbose && enter(verbose, "Non-SNP units:");
      verbose && str(verbose, units);
      if (length(units) > 0) {
        nonSNPs <- getCellIndices(cdf, units=units, 
                       useNames=FALSE, unlist=TRUE, verbose=verbose);
      } else {
        nonSNPs <- NULL;
      }
      # Not needed anymore
      units <- NULL;
      verbose && enter(verbose, "Non-SNP cells:");
      verbose && str(verbose, nonSNPs);
      verbose && exit(verbose);

      cells <- getAlleleCellPairs(cdf, verbose=verbose);
      S <- as.integer(B %/% 2 + 1);
      shifts <- -(S-1):(S-1);
      verbose && cat(verbose, "Probe shifts:");
      verbose && print(verbose, shifts);
      snps <- groupBySnpNucleotides(acs, cells=cells, shifts=shifts, 
                                                     verbose=verbose);
      # Not needed anymore
      cells <- shifts <- NULL;
      # Clean out empty sets
      for (kk in seq_along(snps)) {
        cells <- snps[[kk]];
        if (length(cells) == 0)
          snps[[kk]] <- NULL;
      }
      setsOfProbes <- list(snps=snps, nonSNPs=nonSNPs);
      # Not needed anymore
      snps <- nonSNPs <- NULL;
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Group by number of nucleotides (B) around SNP position?
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (B == 0) {
      verbose && enter(verbose, "Merging SNP groups");
      snpsT <- list(all=NULL);
      snps <- setsOfProbes$snps;
      for (kk in seq_along(snps)) {
        snpsT$all <- cbind(snpsT$all, snps[[kk]]);
      }
      # Not needed anymore
      snps <- NULL;
      setsOfProbes$snps <- snpsT;
      # Not needed anymore
      snpsT <- NULL;
      verbose && exit(verbose);
    } else if (B == 1) {
      # Nothing to do, default
    } else {
      if (version != 4)
        throw("Not supported.");
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Merge shifts?
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (mergeShifts) {
      # Nothing to do.
    } else {
      # Split by shift
      verbose && enter(verbose, "Grouping probe pairs by their shift relative to SNP");

      # For each set of cells...
      for (gg in seq_along(setsOfProbes)) {
        cells <- setsOfProbes[[gg]];
        groupTag <- names(setsOfProbes)[gg];
        verbose && enter(verbose, sprintf("Set #%d ('%s') of %d", gg, groupTag, length(setsOfProbes)));

        dim <- dim(cells);
        if (is.null(dim) || dim[2] != 2) {
          # Not a set of SNP cell indexes.
          next;
        }

        # Identify the location of the SNP position for each probe pair
        naValue <- as.integer(NA);
        snpPositions <- rep(naValue, length=dim[1]);
        possibleShifts <- as.integer(seq(from=-4, to=+4));  # Hardwired!
        possiblePositions <- 13 + possibleShifts;
        for (pos in possiblePositions) {
          seqsA <- readSequenceMatrix(acs, cells=cells[,1], positions=pos);
          seqsB <- readSequenceMatrix(acs, cells=cells[,2], positions=pos);
          snpPositions[(seqsA != seqsB)] <- pos;
          # Not needed anymore
          seqsA <- seqsB <- NULL;
        }
        shifts <- snpPositions - as.integer(13);

        possibleShifts <- sort(unique(shifts), na.last=TRUE);

        # For each shift...
        subgroups <- vector("list", length(possibleShifts));
        shiftTags <- sprintf("shift=%+d", possibleShifts);
        names(subgroups) <- paste(groupTag, shiftTags, sep=",");

        for (ss in seq_along(possibleShifts)) {
          shift <- possibleShifts[ss];
          if (is.na(shift)) {
            idxs <- which(is.na(shifts));
          } else {
            idxs <- which(shifts == shift);
          }
          cellsSS <- cells[idxs,,drop=FALSE];
          subgroups[[ss]] <- cellsSS;
          # Not needed anymore
          idxs <- cellsSS <- NULL;
        } # for (shift ...)
        # Not needed anymore
        cells <- NULL;

        setsOfProbes[[gg]] <- subgroups;
        # Not needed anymore
        subgroups <- NULL;

        verbose && exit(verbose);
      } # for (gg ...)
      verbose && exit(verbose);
    }
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Fake symmetry by flipping every 2nd (A,B) to (B,A)?
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (fakeSymmetry) {
      verbose && enter(verbose, "Flipping every 2nd (A,B) to (B,A)");
      snps <- setsOfProbes$snps;
      for (kk in seq_along(snps)) {
        idxs <- snps[[kk]];
        cc <- seq(from=1, to=ncol(idxs), by=2);
        idxs[,cc] <- idxs[2:1,cc, drop=FALSE];
        snps[[kk]] <- idxs;
        # Not needed anymore
        cc <- idxs <- NULL;
      } # for (gg ...)
      setsOfProbes$snps <- snps;
      # Not needed anymore
      snps <- NULL;
      verbose && exit(verbose);
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Transpose
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # For each set of cells...
    snps <- setsOfProbes$snps;
    for (kk in seq_along(snps)) {
      snps[[kk]] <- t(snps[[kk]]);
    } # for (gg ...)
    setsOfProbes$snps <- snps;
    # Not needed anymore
    snps <- NULL;

    # Tag the result with a version number
    attr(setsOfProbes, "version") <- version;

    this$.setsOfProbes <- setsOfProbes;

    verbose && exit(verbose);
  }

  setsOfProbes;
}, protected=TRUE) # getSetsOfProbes()



############################################################################
# HISTORY:
# 2008-11-28
# o Now getSetsOfProbes(..., version=0) uses the field '.pairBy' to decide
#   what version to use.
# 2008-09-12
# o Added argument 'fakeSymmetry' to getSetsOfProbes().  Follow up: This
#   seems to work; with 'fakeSymmetry=TRUE' the estimates od the 'origin'
#   are very similar for allele A and allele B.  The estimated backtransform
#   matrix 'Winv' is also quite symmetric.
# 2008-09-06
# o BUG FIX: getSetsOfProbes() for version 4 would return 'nonSNPs' 
#   as units, not cells.  Arghh...
# 2008-09-05
# o BUG FIX: getSetsOfProbes() would return 'nonSNPs' in a sublist.  This
#   caused ACC to skip the offset calibration of non-SNP units.
# 2008-09-02
# o UPDATE: Now getSetsOfProbes() returns a list of two elements 'snps' and
#   'nonSNPs'. The 'snps' element is in turn a list of probe pairs groups.
#   The probe pairs are Kx2 matrices, where the rows are now alleles  A & B.
# 2008-08-31
# o BUG FIX: The allele pairs identified was not correct for GWS arrays.
# o Updated AllelicCrosstalkCalibration to support flavor 'expectile' too.
# 2008-08-30
# o Added argument 'mergeShifts=TRUE' and 'B=1'.  Currently B=0 and B=1 
#   is supported.
# 2008-08-29
# o Added protected getSetsOfProbes().  By overriding this method, other
#   sets of probes can be used.  This function might later also recognize
#   user specified function generating the sets.
# 2008-08-04
# o Added support to fit the genotype "cone" using the 'expectile' package
#   instead of the 'sfit' package. This is controlled by the 'flavor' 
#   argument of the constructor.
# 2008-07-14
# o Now explicitly using matrix(..., byrow=FALSE).
# 2008-05-30
# o BUG FIX: The constructor of AllelicCrosstalkCalibration used 
#   non-defined variable 'verbose'.
# 2008-02-21
# o Now SNPs and CN probes are infered from getUnitTypes(cdf) and no longer
#   from the unit names.
# 2008-02-14
# o Now 'verbose' is passed as a logical argument to fitGenotypeCone().
# 2007-12-01
# o MEMORY OPTIMIZATION: Added clearCache() to AllelicCrosstalkCalibration.
# o BUG FIX: The AllelicCrosstalkCalibration introduced in previous version
#   was broken for 10K (maybe 100K and 500K as well).
# o Now 'subsetToAvg' of AllelicCrosstalkCalibration accepts '-XY' (and
#   '-X' and '-Y') for automatic look up of all units and exclude those 
#   that are on ChrX and ChrY.  Note, '-XY' will work on all chip types,
#   also older ones for which there are no ChrY units.
# o The new constructor argument 'rescaleBy' now sets a "subtag", e.g.
#   'ACC,ra' where the 'ra' indicates that 'rescaleBy=all' was used.
# 2007-11-27
# o Added specific getAsteriskTag() for AllelicCrosstalkCalibration.
# o Added the option to rescale towards a target average of *all* probes.
# o Now rescale() throws an error if length(params$targetAvg) is not 1 or 2.
#   This should never happend, but was added just in case.
# 2007-09-14
# o BUG FIX: Now the target average of non-SNP probes is half of the target
#   average of alleles.
# 2007-09-09
# o Added alpha version of plotBasepair() to AllelicCrosstalkCalibration.
# 2007-09-08
# o Now AllelicCrosstalkCalibration corrects also non-SNP PM cells by 
#   substracting a global offset and rescaling towards target average.
#   The global offset is calculated as the weighted average of all
#   allelic offsets.  This is the simplest way to incorporate a calibration
#   for non-SNP cells.  A more advanced version would be to stratify by
#   middle nucleotide, and use the corresponding estimates from the SNP
#   cells, but that would require information about the middle nucleotide
#   for all cells, which got some overhead.  Hopefully the simpler version
#   is good enough.
# 2007-09-05
# o Now the rescaling can be done either on (yA,yB) separately or on 
#   y=yA+yB.  If targetAvg has two values the former, otherwise the latter.
# o Now AllelicCrosstalkCalibration recognizes argument 'subsetToAvg'.
# o Now process() stores the crosstalk settings and estimated parameters
#   to file. May be useful if one wants to go back and look at the details.
#   One day we might get around to store this information in the CEL file
#   headers.
# o Now process() first fits the crosstalk model for all basepairs, then
#   backtransform the signals, then optional rescale signals to target
#   average, then saves the calibrated signals.
# o SPEED UP: Now getAlleleProbePairs() is only called if data needs to be
#   calibrated, i.e. if already calibrated it is not loaded.
# o CLEAN UP: Now the code of calibrateAllelicCrosstalk() is included here.
# 2007-03-29
# o Now 'targetAvg' defaults to 2200 so that allele A and allele B signals
#   are rescaled to be one the same scale.  If so,  it does not make sense
#   to do background correction afterwards.
# o Added getParameters().
# o Added support for arguments 'targetAvg', 'alpha', 'q', and 'Q'.
# 2006-12-08
# o Now this class inherits from the ProbePreprocessing class.
# o Now this pre-processor output results to probeData/.
# o Renamed from AllelicCrosstalkCalibrator.
# 2006-11-18
# o Removed version and subversion tags, and related functions. 
#   Now getTags() returns the tags of the input data set plus any tags 
#   of this instance.
# 2006-11-02
# o Created from QuantileNormalizer.R.
############################################################################
