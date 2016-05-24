###########################################################################/**
# @set "class=AffymetrixCdfFile"
# @RdocMethod computeAffinities
#
# @title "Calculates probe affinities from sequence"
#
# \description{
#  @get "title".
#
# Adapted from @see "gcrma::compute.affinities" in the \pkg{gcrma} package.
# Attempts to find the tab-separated probe sequence file associated with
# a particular CDF, and matches sequence to probe index in order to assign
# an affinity to each probe.
# }
#
# @synopsis
#
# \arguments{
#   \item{safe}{A @logical argument passed to \code{getProbeSequenceData()}.}
#   \item{force}{If @FALSE, cached results is returned, if available.}
#   \item{verbose}{See @see "R.utils::Verbose".}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @numeric @vector of (log2) probe affinities, of length equal
#  to the total number of features on the array.
# }
#
# @author "KS, HB, MR"
#*/###########################################################################
setMethodS3("computeAffinities", "AffymetrixCdfFile", function(this, safe=TRUE, force=FALSE, verbose=FALSE, ...) {
  isPackageInstalled("gcrma") || throw("Package not installed: gcrma");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # getSeqMatrix() corresponds to calling the following native function
  #   .Call("gcrma_getSeq", seq, PACKAGE="gcrma");
  # However, we do want to avoid to load 'gcrma', because it loads Biostrings
  # which loads IRanges, which cause problems, e.g. trim(). /HB 2009-05-09
  getSeqMatrix <- function(seq, ...) {
    seq <- charToRaw(seq);
    map <- charToRaw("ACGT");
    res <- matrix(0L, nrow=4, ncol=length(seq));
    for (kk in 1:4) {
      res[kk,(seq == map[kk])] <- 1L;
    }
    res;
  } # getSeqMatrix()

  getGcrmaSplineCoefs <- function(...) {
    env <- new.env();
    data("affinity.spline.coefs", package="gcrma", envir=env);
    env$affinity.spline.coefs;
  } # getGcrmaSplineCoefs()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'safe':
  safe <- Arguments$getLogical(safe);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  chipTypeFull <- getChipType(this, fullname=TRUE);


  verbose && enter(verbose, "Computing GCRMA probe affinities");
  verbose && cat(verbose, "Number of units: ", nbrOfUnits(this));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Looking for MMs (and PMs) in the CDF
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identify PMs and MMs among the CDF cell indices");
  isPm <- isPm(this);
  verbose && str(verbose, isPm);
  verbose && summary(verbose, isPm);
  nbrOfPMs <- sum(isPm);
  verbose && cat(verbose, "MMs are defined as non-PMs");
  isMm <- !isPm;
  nbrOfMMs <- sum(isMm);
  verbose && cat(verbose, "Number of PMs: ", nbrOfPMs);
  verbose && cat(verbose, "Number of MMs: ", nbrOfMMs);
  verbose && exit(verbose);

  # Sanity check
  if (nbrOfMMs == 0) {
#    throw("Cannot calculate gcRMA probe affinities. The CDF contains no MMs: ", chipTypeFull);
  }

  # Checking cache
  key <- list(method="computeAffinities", class=class(this)[1],
              chipTypeFull=chipTypeFull, version="2009-05-09");
  if (getOption(aromaSettings, "devel/useCacheKeyInterface", FALSE)) {
    key <- getCacheKey(this, method="computeAffinities", chipTypeFull=chipTypeFull);
  }
  dirs <- c("aroma.affymetrix", chipTypeFull);
  if (!force) {
    res <- loadCache(key=key, dirs=dirs);
    if (!is.null(res)) {
      verbose && cat(verbose, "Cached results found.");
      verbose && exit(verbose);
      return(res);
    }
  }


  verbose && enter(verbose, "Reading probe-sequence data");
  sequenceInfo <- getProbeSequenceData(this, safe=safe, verbose=verbose);
  verbose && str(verbose, sequenceInfo);
  # Dropping unused (x,y) columns (not really needed, but towards using ACS files)
  sequenceInfo <- sequenceInfo[,c("cell", "sequence"), drop=FALSE];
  # Sorting by cell index (not really needed, but towards using ACS files)
  o <- order(sequenceInfo$cell);
  sequenceInfo <- sequenceInfo[o,,drop=FALSE];
  verbose && str(verbose, sequenceInfo);
  nbrOfSequences <- nrow(sequenceInfo);
  verbose && printf(verbose, "Number of known sequences: %d of %d (%.1f%%)\n",
               nbrOfSequences, nbrOfCells(this), 100*nbrOfSequences/nbrOfCells(this));
  verbose && exit(verbose);


  verbose && enter(verbose, "Getting all CDF cell indices for all units");
  cells <- getCellIndices(this, useNames=FALSE, unlist=TRUE, verbose=verbose);
  verbose && str(verbose, cells);
  n <- length(cells);
  verbose && printf(verbose, "Number of CDF cell indices: %d of %d (%.1f%%)\n",
                                   n, nbrOfCells(this), 100*n/nbrOfCells(this));
  verbose && exit(verbose);

  sum(is.element(sequenceInfo$cell, cells));
  verbose && printf(verbose, "Number of known sequences out of all CDF cell indices: %d of %d (%.1f%%)\n",
                             nbrOfSequences, n, 100*nbrOfSequences/n);


  # probe sequence tab-separated files generally only contain the PM
  # sequence (since we know that MM sequence always has base 13 complemented).
  # This is true for expression arrays, but possibly not for genotyping arrays.
  # We will calculate affinity for PM and MM anyway, and then
  # try to match MM indices from the CDF file with their appropriate PMs (and
  # hence assign the correct affinity).

  # assume that MM has same x-coordinate as PM, but y-coordinate larger by 1
  indexPm <- sequenceInfo$cell;
  dimension <- getDimension(this);
  indexMmPutative <- indexPm + dimension[1];

  # match up putative MM with actual list of MMs from CDF
  matches <- match(indexMmPutative, cells[isMm]);
  indexMm <- cells[isMm][matches];
  notNA <- which(!is.na(indexMm));
  indexMm <- indexMm[notNA];

  # for PM+MM arrays, the number of MMs and the number of PMs in the
  # CDF is equal
  isPMMMChip <- (length(indexMm) == length(indexPm));
  verbose && enter(verbose, "Chip type is an \"PM+MM\" chip: ", isPMMMChip);

  # Not needed anymore
  indexMmPutative <- matches <- NULL; # Not needed anymore


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate basis matrix based on priori probe affinitites (of gcrma)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # now calculate affinities - code reused from compute.affinities() in gcrma
  verbose && enter(verbose, "Calculating probe affinities");

  # To please R CMD check on R v2.6.0
  affinity.spline.coefs <- getGcrmaSplineCoefs();

  affinity.basis.matrix <- splines::ns(1:25, df=length(affinity.spline.coefs)/3);

  A13 <- sum(affinity.basis.matrix[13, ] * affinity.spline.coefs[1:5]);
  T13 <- 0;
  C13 <- sum(affinity.basis.matrix[13, ] * affinity.spline.coefs[6:10]);
  G13 <- sum(affinity.basis.matrix[13, ] * affinity.spline.coefs[11:15]);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Predicting probe affinities based on the probe sequences
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfCells <- dimension[1]*dimension[2];
  # Not needed anymore
  dimension <- NULL; # Not needed anymore

  # Allocate empty vector of affinities
  naValue <- as.double(NA);
  affinities <- rep(naValue, times=nbrOfCells);

  if(isPMMMChip) {
    # ---------------------------------------------------------
    # keep this code below untouched for the PM+MM expression arrays
    # ---------------------------------------------------------
    apm <- vector("double", nbrOfSequences);
    amm <- vector("double", nbrOfSequences);

    # Differencies in affinities between the PM and the MM probe
    deltas <- c(A=T13-A13, C=G13-C13, G=C13-G13, T=A13-T13);

    sequences <- sequenceInfo$sequence;
    # Not needed anymore
    sequenceInfo <- NULL;

    pb <- NULL;
    if (verbose) {
      cat(verbose, "Progress (counting to 100): ");
      if (isVisible(verbose)) {
        pb <- ProgressBar(stepLength=100/(nbrOfSequences/1000));
        reset(pb);
      }
    }

    for (ii in seq_along(apm)) {
      if (!is.null(pb) && (ii %% 1000 == 0)) {
        increase(pb);
      }

      # Get a 4x25 matrix with rows A, C, G, and T.
      charMtrx <- getSeqMatrix(sequences[ii]);

      # Calculate the PM affinity
      A <- cbind(charMtrx[1,,drop=TRUE] %*% affinity.basis.matrix,
                 charMtrx[2,,drop=TRUE] %*% affinity.basis.matrix,
                 charMtrx[3,,drop=TRUE] %*% affinity.basis.matrix);
      apm[ii] <- A %*% affinity.spline.coefs;

      # Calculate the MM affinity from the PM affinity by
      # correcting the mismatch nucleotide.
      base <- which(charMtrx[,13,drop=TRUE] == 1);
      # Sanity check
      if (length(beta) != 1) throw("Unknown nucleotide index: ", base);
      amm[ii] <- apm[ii] + deltas[base];
    } # for (ii in ...)
    # Not needed anymore
    charMtrx <- A <- NULL; # Not needed anymore

    if (!is.null(pb)) {
      increase(pb);
    }

    # create a vector to hold affinities and assign values to the
    # appropriate location in the vector
    affinities[indexPm] <- apm;
    affinities[indexMm] <- amm[notNA];
    # Not needed anymore
    indexPm <- indexMm <- apm <- amm <- notNA <- NULL; # Not needed anymore
    # ---------------------------------------------------------
  } else {
    # ---------------------------------------------------------
    # new code to compute affinities from the MM probes
    # (antigenomic or whatever) on PM-only arrays  -- MR 2009-03-28
    # ---------------------------------------------------------
    naValue <- as.double(NA);
    apm <- rep(naValue, times=nbrOfSequences);

    indexAll <- sequenceInfo$cell;
    sequences <- sequenceInfo$sequence;
    # Not needed anymore
    sequenceInfo <- NULL;

    # Process only sequences of length 25
    n <- nchar(sequences);
    idxs <- which(n == 25);
    nbrOfNon25mers <- nbrOfSequences - length(idxs);
    if (nbrOfNon25mers > 0) {
      apm[-idxs] <- naValue;  # It is already NA?!? /HB 2010-09-29
      warning("Detected ", nbrOfNon25mers, " sequence that are not of length 25 nucleotides. For those probes, the affinities are defined to be NA.");
    }

    pb <- NULL;
    if (verbose) {
      cat(verbose, "Progress (counting to 100): ");
      if (isVisible(verbose)) {
        pb <- ProgressBar(stepLength=100/(length(idxs)/1000));
        increase(pb);
      }
    }

    for (ii in seq_along(idxs)) {
      if (!is.null(pb) && (ii %% 1000 == 0)) {
        increase(pb);
      }
      idx <- idxs[ii];
      # Get a 4x25 matrix with rows A, C, G, and T.
      charMtrx <- getSeqMatrix(sequences[idx]);

      # Calculate the PM affinity
      A <- cbind(charMtrx[1,,drop=TRUE] %*% affinity.basis.matrix,
                 charMtrx[2,,drop=TRUE] %*% affinity.basis.matrix,
                 charMtrx[3,,drop=TRUE] %*% affinity.basis.matrix);
      apm[idx] <- A %*% affinity.spline.coefs;
    } # for (ii in ...)

    if (!is.null(pb)) {
      increase(pb);
    }

    affinities[indexAll] <- apm;

    # Not needed anymore
    indexAll <- apm <- NULL; # Not needed anymore
  } # if(isPMMMChip)
  verbose && exit(verbose);

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  # Sanity checks
  stopifnot(length(affinities) == nbrOfCells(this));

  # Saving to cache
  comment <- paste(unlist(key, use.names=FALSE), collapse=";");
  saveCache(key=key, affinities, comment=comment, dirs=dirs);

  verbose && exit(verbose);

  affinities;
}, private=TRUE)



setMethodS3("getACSFile", "AffymetrixCdfFile", function(this, ...) {
  getAromaCellSequenceFile(this, ...);
})


setMethodS3("getAromaCellSequenceFile", "AffymetrixCdfFile", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Locating the aroma cell sequence (ACS) annotation data file");
  chipTypeS <- getChipType(this, fullname=FALSE);
  verbose && cat(verbose, "Chip type: ", chipTypeS);
  nbrOfCells <- nbrOfCells(this);
  verbose && cat(verbose, "Number of cells: ", nbrOfCells);
  acs <- AromaCellSequenceFile$byChipType(chipTypeS, nbrOfCells=nbrOfCells);
  verbose && print(verbose, acs);
  verbose && exit(verbose);

  acs;
})


setMethodS3("computeAffinitiesByACS", "AffymetrixCdfFile", function(this, ..., method=c("v3", "v2", "v1"), force=FALSE, verbose=FALSE) {
  isPackageInstalled("gcrma") || throw("Package not installed: gcrma");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # getSeqMatrix() corresponds to calling the following native function
  #   .Call("gcrma_getSeq", seq, PACKAGE="gcrma");
  # However, we do want to avoid to load 'gcrma', because it loads Biostrings
  # which loads IRanges, which cause problems, e.g. trim(). /HB 2009-05-09
  getSeqMatrix <- function(seq, ...) {
    seq <- charToRaw(seq);
    map <- charToRaw("ACGT");
    res <- matrix(0L, nrow=4, ncol=length(seq));
    for (kk in 1:4) {
      res[kk,(seq == map[kk])] <- 1L;
    }
    res;
  } # getSeqMatrix()

  getGcrmaSplineCoefs <- function(...) {
    env <- new.env();
    data("affinity.spline.coefs", package="gcrma", envir=env);
    env$affinity.spline.coefs;
  } # getGcrmaSplineCoefs()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'method':
  method <- match.arg(method);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  verbose && enter(verbose, "Computing GCRMA probe affinities");
  chipTypeS <- getChipType(this, fullname=FALSE);
  verbose && cat(verbose, "Chip type: ", chipTypeS);
  verbose && cat(verbose, "Number of units: ", nbrOfUnits(this));


  # Checking cache
  key <- list(method="computeAffinities", class=class(this)[1],
              chipType=chipTypeS, version="2010-09-29");
  if (getOption(aromaSettings, "devel/useCacheKeyInterface", FALSE)) {
    key <- getCacheKey(this, method="computeAffinitiesByACS", chipType=chipTypeS);
  }
  dirs <- c("aroma.affymetrix", chipTypeS);
  if (!force) {
    res <- loadCache(key=key, dirs=dirs);
    if (!is.null(res)) {
      verbose && cat(verbose, "Cached results found.");
      verbose && exit(verbose);
      return(res);
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Locate the cell sequence annotation data file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Locating the cell sequence annotation data file");
  acs <- getACSFile(this);
  verbose && print(verbose, acs);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate basis matrix based on priori probe affinitites (of gcrma)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Calculating apriori gcRMA affinities parameters");
  # now calculate affinities - code reused from compute.affinities() in gcrma
  verbose && enter(verbose, "Calculating probe affinities");

  # To please R CMD check on R v2.6.0
  affinity.spline.coefs <- getGcrmaSplineCoefs();
  df <- length(affinity.spline.coefs) / 3;
  affinity.basis.matrix <- splines::ns(1:25, df=df);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identifying cells with known probe sequences
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identifying cells with known probe sequences");
  nbrOfCells <- nbrOfCells(this);
  bp <- readSequenceMatrix(acs, positions=1, what="character");
  hasSequence <- !is.na(bp);
  cells <- which(hasSequence);
  verbose && printf(verbose, "Number of known sequences: %d of %d (%.1f%%)\n",
                     length(cells), nbrOfCells, 100*length(cells)/nbrOfCells);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Predicting probe affinities based on the probe sequences
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Calculating chip-type specific gcRMA probe affinities");

  # Allocate empty vector of affinities
  naValue <- as.double(NA);
  affinities <- rep(naValue, times=nbrOfCells);

  if (method == "v1") {
    pb <- NULL;
    if (verbose) {
      cat(verbose, "Progress (counting to 100): ");
      if (isVisible(verbose)) {
        pb <- ProgressBar(stepLength=100/(length(cells)/1000));
        reset(pb);
      }
    }

    verbose && enter(verbose, "Reading all cell sequences");
    seqMatrix <- readSequenceMatrix(acs, cells=cells);
    verbose && str(verbose, seqMatrix);
    verbose && exit(verbose);

    bases <- c("A", "C", "G", "T");
    for (ii in seq_along(cells)) {
      if (!is.null(pb) && (ii %% 1000 == 0)) {
        increase(pb);
      }

      seqMatrixII <- seqMatrix[ii,,drop=TRUE];
      affinity <- 0;
      for (kk in 1:3) {
        base <- bases[kk];
        idxs <- (df*(kk-1)) + 1:5;
        coefsKK <- affinity.spline.coefs[idxs];
        baseKK <- as.integer(seqMatrixII == base);
        Akk <- baseKK %*% affinity.basis.matrix;
        akk <- sum(Akk * coefsKK);
        affinity <- affinity + akk;
      }

      cell <- cells[ii];
      affinities[cell] <- affinity;
    } # for (ii ...)
    if (!is.null(pb)) {
      increase(pb);
    }
  } else if (method == "v2") {
    pb <- NULL;
    if (verbose) {
      cat(verbose, "Progress (counting to 100): ");
      if (isVisible(verbose)) {
        pb <- ProgressBar(stepLength=100/(length(cells)/1000));
        reset(pb);
      }
    }

    verbose && enter(verbose, "Reading all cell sequences");
    seqMatrix <- readSequenceMatrix(acs, cells=cells);
    verbose && str(verbose, seqMatrix);
    verbose && exit(verbose);

    bases <- c("A", "C", "G", "T");
    affinitiesT <- double(length(cells));
    # For each nucleotide type...
    for (kk in 1:3) {
      base <- bases[kk];
      idxs <- (df*(kk-1)) + (1:df);
      coefsKK <- affinity.spline.coefs[idxs];

      # For each sequence...
      for (ii in seq_along(cells)) {
        seqMatrixII <- seqMatrix[ii,,drop=TRUE];
        isBase <- (seqMatrixII == base);
        basesII <- which(isBase);

        # For each position...
        Bkk <- double(df);
        for (pp in basesII) {
          Bpp <- affinity.basis.matrix[pp,,drop=TRUE];
          Bkk <- Bkk + Bpp;
        }
        b <- sum(Bkk * coefsKK);
        affinitiesT[ii] <- affinitiesT[ii] + b;
      } # for (ii ...)
    } # for (kk ...)
    affinities[cells] <- affinitiesT;
##    print(affinities[cells] == affinities1[cells]);
    if (!is.null(pb)) {
      increase(pb);
    }
  } else if (method == "v3") {
    nbrOfPositions <- getProbeLength(acs);
    bases <- c("A", "C", "G", "T");
    affinitiesT <- double(length(cells));
    # For each position...
    for (pp in seq_len(nbrOfPositions)) {
      verbose && enter(verbose, sprintf("Nucleotide position #%d of %d", pp, nbrOfPositions));

      affinityBasisPP <- affinity.basis.matrix[pp,,drop=TRUE];
      verbose && cat(verbose, "Affinity basis at this position:");
      verbose && print(verbose, affinityBasisPP);

      # Nothing to do?
      if (all(affinityBasisPP == 0)) {
        verbose && exit(verbose);
        next;
      }

      # The nucleotides at this position across all cells
      verbose && enter(verbose, "Reading cell nucleotides at this position");
      nucleotides <- readSequenceMatrix(acs, cells=cells, positions=pp, drop=TRUE);
      verbose && str(verbose, nucleotides);
      verbose && exit(verbose);

      # For each nucleotide type...
      for (tt in 1:3) {
        base <- bases[tt];
        verbose && enter(verbose, "Nucleotide ", base);

        # Identify cells with this nucleotide at this position
        idxsTT <- which(nucleotides == base);

        verbose && cat(verbose, "Indices of cells with this base:");
        verbose && str(verbose, idxsTT);

        # Nothing to do?
        if (length(idxsTT) == 0) {
          verbose && exit(verbose);
          next;
        }

        idxs <- (df*(tt-1)) + (1:df);
        coefsTT <- affinity.spline.coefs[idxs];
        affinityBasisPPTT <- affinityBasisPP * coefsTT;
        affinityPPTT <- sum(affinityBasisPPTT);

        verbose && printf(verbose, "Affinity for nucleotide '%s' at position #%d: %.3f\n", base, pp, affinityPPTT);

        affinitiesT[idxsTT] <- affinitiesT[idxsTT] + affinityPPTT;

        verbose && exit(verbose);
      } # for (tt ...)

      verbose && exit(verbose);
    } # for (pp ...)
    affinities[cells] <- affinitiesT;
    ##  print(isZero(affinities[cells] - affinities1[cells], neps=5));
  } # if (method ...)

  verbose && exit(verbose);

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  # Saving to cache
  comment <- paste(unlist(key, use.names=FALSE), collapse=";");
  saveCache(key=key, affinities, comment=comment, dirs=dirs);

  verbose && exit(verbose);

  affinities;
}, private=TRUE)


############################################################################
# HISTORY:
# 2014-06-29
# o BUG FIX: getAromaCellSequenceFile() for AffymetrixCdfFile used
#   undefined variable 'nbrOfCells'.
# 2010-10-01
# o Added getAromaCellSequenceFile() and getACSFile() for AffymetrixCdfFile.
# 2010-09-30
# o Added computeAffinitiesByACS() for AffymetrixCdfFile, which calculated
#   gcRMA affinities straight off from cell sequences in ACS file.
#   It does not have to do funny PM & MM lookups and it is much faster!
# 2010-09-29
# o Dropped non-used argument 'path' from computeAffinities().
# o Now the progress bar produced by computeAffinities() is ended correctly.
# 2010-04-15
# o BUG FIX: computeAffinities(..., verbose=FALSE) of AffymetrixCdfFile
#   would give throw "Error in reset(pb) : object 'pb' not found".
#   Thanks Stephen ? at Mnemosyne BioSciences, Finland.
# 2009-05-09 [HB]
# o CLEAN UP: computeAffinities() of AffymetrixCdfFile no longer need to
#   load 'gcrma'.  However, it still needs to load a small set of
#   prior estimates from the 'gcrma' package.
#   'matchprobes' package, but never used it.
# o CLEAN UP: computeAffinities() of AffymetrixCdfFile did load the
#   'matchprobes' package, but never used it.
# o Added private getProbeSequenceData() to AffymetrixCdfFile. This method
#   will later retrieve probe sequences from the binary ACS file instead
#   of the probe-tab files.
# o Updated to make full use of the AffymetrixProbeTabFile class, which
#   translateColumnNames() will take care of all variations of names for
#   the column containing unit names.
# 2009-03-28 [MR]
# o Made several modifications to allow computing of affinities for
#   Gene 1.0 ST arrays.  For example:
#   -- left the PM+MM array code mostly untouched
#   -- fixed some assumptions about the columns of the probe_tab file
#   -- added a different stream for PM-only (with NCs)
# 2007-07-30
# o UPDATE: Now computeAffinities() for AffymetrixCdfFile gives an error
#   if there are no MMs in the CDF.
# o BUG FIX: computeAffinities() for AffymetrixCdfFile searched for the
#   probe-tab file using the chip type given by the fullname of the CDF
#   and not the basic name.
# 2007-09-06
# o Made computeAffinities() more memory efficient since it is using the
#   new getCellIndices(cdf, useNames=FALSE, unlist=TRUE).
# 2007-06-07
# o BUG FIX: If an Affymetrix probe tab file is not found for the chip type,
#   computeAffinitities() of AffymetrixCdfFile would throw "Error in
#   paste(..., sep = sep, collapse = collapse): object "pattern" not found"
#   instead of an intended and more information error.
# 2007-02-06
# o Now computeAffinities() is locating the Affymetrix' probe-sequence file
#   using AffymetrixProbeTabFile$findByChipType().
# 2006-10-09
# o Added caching to file.
# o Now default probe affinity is NA (for non estimated probes).
# o Added a progress bar to the calculations.
# o Made internal reader smart so it guess the X, Y, and sequence columns
#   by reading the data and comparing to CDF information.
# o Added verbose statements.
# o Made sure connection is closed regardless how the method exits.
# 2006-10-04
# o Debugged, tested, docs/comments added.
# 2006-10-01
# o Created.
############################################################################
