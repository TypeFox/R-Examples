###########################################################################/**
# @RdocClass MatSmoothing
#
# @title "The MatSmoothing class"
#
# \description{
#  @classhierarchy
#
#  This class represents a function for smoothing data with a trimmed mean.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "ProbeLevelTransform".}
#   \item{design}{A design @matrix.}
#   \item{probeWindow}{Bandwidth to use.  Effectively the width is
#      2*probeWindow since it looks probeWindow bases in either direction.}
#   \item{nProbes}{The minimum number of probes to calculate a MAT score for.}
#   \item{meanTrim}{The amount of trimming of the mean in [0,0.5].}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "MR, HB"
#*/###########################################################################
setConstructorS3("MatSmoothing", function(..., design=NULL, probeWindow=300, nProbes=10, meanTrim=0.1) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(design)) {
    # Argument 'design':
    if (!is.matrix(design)) {
      throw("Argument 'design' is not a matrix: ", class(design)[1]);
    }

    # Argument 'probeWindow':
    probeWindow <- Arguments$getNumeric(probeWindow, range=c(0,Inf));

    # Argument 'nProbes':
    nProbes <- Arguments$getInteger(nProbes, range=c(1,Inf));

    # Argument 'meanTrim':
    meanTrim <- Arguments$getNumeric(meanTrim, range=c(0,0.5));
  }


  this <- extend(ProbeLevelTransform(...), "MatSmoothing",
    .design = design,
    .probeWindow = probeWindow,
    .nProbes = nProbes,
    .meanTrim = meanTrim
  );


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Further validation of the design matrix
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ds <- getInputDataSet(this);
  if (!is.null(ds)) {
    # Validate the dimension of the design matrix
    nbrOfFiles <- length(ds);
    design <- this$.design;
    dim <- dim(design);
    if (dim[1] != nbrOfFiles) {
      throw("The number of rows in the 'design' matrix, does not match the number of arrays in the input data set: ", dim[1], " != ", nbrOfFiles);
    }

#    if (dim[1] != dim[2]) {
#      throw("The 'design' matrix must be a square matrix: ", dim[1], "x", dim[2]);
#    }

    # Validate the contents of the design matrix
    for (cc in seq_len(ncol(design))) {
      if (!any(design[,cc] != 0)) {
        throw("Column #", cc, " in argument 'design' is all zero.");
      }
    } # for (cc ...)

    # Validate the names attributes of the design matrix
    outputNames <- colnames(design);
    if (is.null(outputNames)) {
      throw("Matrix 'design' does not have column names.");
    }
    outputNames <- Arguments$getCharacters(outputNames);
    if (any(duplicated(outputNames))) {
      throw("Argument 'design' contains duplicated column names: ",
             paste(outputNames[duplicated(outputNames)]), collapse=", ");
    }

    # Check if the column names translates to valid filenames
    fullnames <- lapply(outputNames, FUN=function(fullname) {
      Arguments$getFilename(fullname);
    });
  }

  this;
})


setMethodS3("getAromaCellPositionFile", "MatSmoothing", function(this, ..., force=FALSE) {
  acp <- this$.acp;

  if (force || is.null(acp)) {
    dataSet <- getInputDataSet(this);
    cdf <- getCdf(dataSet);
    chipType <- getChipType(cdf, fullname=FALSE);
    nbrOfCells <- nbrOfCells(cdf);
    acp <- AromaCellPositionFile$byChipType(chipType, nbrOfCells=nbrOfCells, ...);
    this$.acp <- acp;
  }

  acp;
}, protected=TRUE)


setMethodS3("getParameters", "MatSmoothing", function(this, ...) {
  # Get parameters from super class
  params <- NextMethod("getParameters");

  # Get parameters of this class
  params2 <- list(
    design = this$.design,
    probeWindow = this$.probeWindow,
    nProbes = this$.nProbes,
    meanTrim = this$.meanTrim
  );

  # Append the two sets
  params <- c(params, params2);

  params;
}, protected=TRUE)



setMethodS3("getExpectedOutputFullnames", "MatSmoothing", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Match to the column names of the design matrix");
  params <- getParameters(this);
  design <- params$design;

  verbose && cat(verbose, "Expected result names:");
  fullnames <- colnames(design);
  verbose && str(verbose, fullnames);

  # Sanity check (backup)
  stopifnot(!is.null(fullnames));

  verbose && exit(verbose);

  fullnames;
}, protected=TRUE)




###########################################################################/**
# @RdocMethod process
#
# @title "Processes the data set"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{force}{If @TRUE, data already processed is re-processed,
#       otherwise not.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @double @vector.
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("process", "MatSmoothing", function(this, ..., units=NULL, force=FALSE, verbose=FALSE) {
  requireNamespace("gsmoothr") || throw("Package not loaded: gsmoothr")
  tmeanC <- gsmoothr::tmeanC;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # internal function to do trimmed mean smoothing
  # Lookup gsmoothr::tmeanC() once; '::' is expensive
  gsmoothr_tmeanC <- gsmoothr::tmeanC;
  calcSmoothed <- function(posVector, dataMatrix, probeWindow, nProbes, meanTrim) {
    nc <- ncol(dataMatrix);
    posM <- matrix(posVector, nrow=length(posVector), ncol=nc);
    o <- order(posM);  # calculate ordering

    smoothedScore <- gsmoothr_tmeanC(posM[o], dataMatrix[o],
           probeWindow=probeWindow, nProbes=nProbes*nc, trim=meanTrim);
    subsetInd <- seq(from=1, to=length(o), by=nc);
    return(smoothedScore[subsetInd]);
  } # calcSmoothed()


  # calculates null distribution from first taking non-overlapping probes,
  # then replicating the negative half of the distribution.
  calcNullDist0 <- function(ch, ps, x) {
    MIN <- -999999
    #inds <-
    y <- rep(MIN, times=length(x))
    n <- length(ch)
    indices <- split(seq_len(n), ch)
    nChr <- length(indices)
    count <- 0
    for (ii in seq_len(nChr)) {
      ind <- indices[[ii]]
      nInd <- length(ind)
      pos <- ind[1]
      for (jj in seq_len(nInd)) {
        pp <- ps[ ind[jj] ]
        if ( (pp-pos) > (probeWindow*2) ) {
          count <- count + 1
          y[count] <- x[ ind[jj] ]
          #inds[count] <- ind[jj]
          pos <- pp
        }
      } # for (jj ...)
    } # for (ii ...)
    y <- y[y > MIN]
    md <- median(y)
    y <- y[y <= md]
    #list( m=md, sd = sd( c(y,-y+2*md) ), inds=inds[inds > MIN])
    list( m=md, sd = sd( c(y,-y+2*md) ) )
  } # calcNullDist0()

  # A version that is almost twice as fast. /HB 2009-05-27
  calcNullDist <- function(chr, pos, x) {
    dblProbeWindow <- 2*probeWindow;

    n <- length(chr);
    indices <- split(seq_len(n), chr);
    nbrOfChromosomes <- length(indices);

    # For each chromosome
    yList <- list();
    for (ii in seq_len(nbrOfChromosomes)) {
      ind <- indices[[ii]];
      nInd <- length(ind);

      # Extract positions
      posII <- pos[ind];
      xII <- x[ind];

      # Reorder along chromosome (don't assume)
      o <- order(posII);
      posII <- posII[o];
      xII <- xII[o];

      yII <- rep(-Inf, times=nInd);

      # For each position on the current chromosome
      lastPos <- posII[1];
      for (jj in seq_len(nInd)) {
        # Find next position outside the probe window
        posJJ <- posII[jj];
        if (posJJ-lastPos > dblProbeWindow) {
          yII[jj] <- xII[jj];
          lastPos <- posJJ;
        }
      } # for (jj ...)

      yList[[ii]] <- yII;
    } # for (ii ...)
    y <- unlist(yList, use.names=FALSE);

    y <- y[y > -Inf];
    md <- median(y);
    y <- y[y <= md];
    list(m=md, sd=sd(c(y,-y+2*md)));
  } # calcNullDist()



  # A version that is even faster if ncol(X) > 1.
  calcNullDists <- function(chr, pos, X) {
    dblProbeWindow <- 2*probeWindow;

    K <- ncol(X);
    names <- colnames(X);

    n <- length(chr);
    idxList <- split(seq_len(n), chr);

    # For each chromosome
    YList <- list();
    for (ii in seq_along(idxList)) {
      idxs <- idxList[[ii]];

      # Extract ordered along chromosome (don't assume)
      posII <- pos[idxs];
      o <- order(posII);
      idxs <- idxs[o];

      posII <- posII[o];
      XII <- X[idxs,,drop=FALSE];

      J <- length(idxs);
      YII <- matrix(-Inf, nrow=J, ncol=K);

      # For each position on the current chromosome
      lastPos <- posII[1];
      for (jj in seq_len(J)) {
        # Find next position outside the probe window
        posJJ <- posII[jj];
        if (posJJ-lastPos > dblProbeWindow) {
          YII[jj,] <- XII[jj,];
          lastPos <- posJJ;
        }
      } # for (jj ...)

      YList[[ii]] <- YII;
    } # for (ii ...)

    m <- sd <- double(K);
    names(m) <- names(sd) <- names;

    res <- list();
    for (kk in seq_len(K)) {
      yList <- lapply(YList, FUN=function(Y) Y[,kk, drop=TRUE]);
      y <- unlist(yList, use.names=FALSE);
      y <- y[y > -Inf];
      m[kk] <- median(y, na.rm=TRUE);
      y <- y[y <= m[kk]];
      sd[kk] <- sd(c(y,-y+2*m[kk]), na.rm=TRUE);
    }

    list(m=m, sd=sd);
  } # calcNullDists()



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'units':
  if (!is.null(units)) {
    units <- Arguments$getIndices(units, max=nbrOfUnits(cdf));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "MAT smoothing according to design matrix");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Already done?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!force && isDone(this)) {
    verbose && cat(verbose, "Already smoothed");
    verbose && exit(verbose);
    outputDataSet <- getOutputDataSet(this);
    return(invisible(outputDataSet));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get input data set
  ds <- getInputDataSet(this);

  # Get (and create) the output path
  outputPath <- getPath(this);

  # Get cdf
  cdf <- getCdf(ds);

  # Get which units to fit
  if (is.null(units)) {
    units <- seq_len(nbrOfUnits(cdf));
  }

  verbose && enter(verbose, "Locating probe position data");
  # Locate AromaCellPositionFile holding probe sequences
  acp <- getAromaCellPositionFile(this, verbose=less(verbose, 5));
  verbose && exit(verbose);

  # Get algorithm parameters
  params <- getParameters(this, verbose=less(verbose, 50));
  probeWindow <- params$probeWindow;
  design <- params$design;
  nProbes <- params$nProbes;
  meanTrim <- params$meanTrim;
  # Not needed anymore
  params <- NULL;


  # ------------------------------------------------------
  # get CEL indices for units
  # ------------------------------------------------------
  # Sanity check to validate that each group for tiling array
  # has this 1 element
  nbrGroupsPerUnit <- nbrOfGroupsPerUnit(cdf);
  stopifnot(all(nbrGroupsPerUnit == 1));

  verbose && enter(verbose, "Loading cell PM indices structured according to the CDF");
  cdfCellsList <- getCellIndices(cdf, units=units, stratifyBy="pm");
  unitNames <- names(cdfCellsList);
  names(cdfCellsList) <- NULL;
  verbose && exit(verbose);

  cellsList <- lapply(cdfCellsList, FUN=function(unit) {
    unit$groups[[1]]$indices;
  });

  nRows <- sapply(cellsList, FUN=length);
  allInds <- unlist(cellsList, use.names=FALSE);

  nbrOfUnits <- length(units);

  fullnamesOut <- colnames(design);
  # Sanity check (backup)
  stopifnot(!is.null(fullnamesOut));
  verbose && cat(verbose, "Result/output names:");
  verbose && str(verbose, fullnamesOut);


  # Preload all cell positions in order to avoid querying
  # the annotation file for each unit. /HB 2009-05-27
  verbose && enter(verbose, "Preloading genomic positions for all cells");
  acpData <- acp[,1:2,drop=FALSE];
  verbose && exit(verbose);

  # --------------------------------------------------------
  # loop through the number of columns in design matrix
  # calculate smoothed score for each combination of samples
  # --------------------------------------------------------
  for (ii in seq_len(ncol(design))) {
    fullname <- fullnamesOut[ii];
    verbose && enter(verbose, sprintf("Result file #%d ('%s') of %d",
                                                  ii, fullname, ncol(design)));

    filename <- sprintf("%s.CEL", fullname);
    pathname <- Arguments$getWritablePathname(filename, path=outputPath, ...);

    # Already done?
    if (!force && isFile(pathname)) {
      verbose && cat(verbose, "Already processed.");
      verbose && exit(verbose);
      next;
    }

    # Check already here if there is already a temporary output file.  This
    # may be a left over from a interrupted previous run, or the fact that
    # the array is processed by another session elsewhere. /HB
    # (allow rename of existing one if forced)
    isFile <- (force && isFile(pathname));
    pathnameT <- pushTemporaryFile(pathname, isFile=isFile, verbose=verbose);

    matScoreNeg <- matScorePos <- outputList <- vector("list", nbrOfUnits);
    names(outputList) <- unitNames;

    sampsKeep <- which(design[,ii] != 0);

    # Sanity check (already done in the setup checks, but a 2nd backup)
    stopifnot(length(sampsKeep) > 0);

    posOrNeg <- design[sampsKeep,ii];

    verbose && enter(verbose, "Reading probe data for the samples needed");
    dsII <- extract(ds, sampsKeep);
    # Reuse the above cdfIndices structure...
    dataList <- readUnits(dsII, units=cdfCellsList, verbose=verbose);
    # ...instead of rereading all again
    ## dataList <- readUnits(dsII, units=units, stratifyBy="pm", verbose=verbose);
    # Not needed anymore
    dsII <- NULL;
    verbose && exit(verbose);

    dataList <- lapply(dataList, FUN=function(u) {
      matrix(log2(u[[1]]$intensities), ncol=length(sampsKeep))
    });

    verbose && enter(verbose, "Computing trimmed means for all units");

    # ------------------------------------------------------
    # loop through each unit
    # ------------------------------------------------------
    for (jj in seq_len(nbrOfUnits)) {
      # allocate space first time through
      if (is.null(outputList[[jj]])) {
        zeroes <- rep(0, times=nRows[jj]);
        matScorePos[[jj]] <- zeroes;
        matScoreNeg[[jj]] <- zeroes;
        outputList[[jj]] <- zeroes;
      }
      if (jj %% 1000 == 0) {
        verbose && cat(verbose, sprintf("Completed %d out of %d units...", jj, nbrOfUnits));
      }

      # loop through all columns in matrix
      sampPos <- which(posOrNeg > 0);
      sampNeg <- which(posOrNeg < 0);
      nPos <- length(sampPos);
      nNeg <- length(sampNeg);

      # extract probe positions
      cells <- cellsList[[jj]];
      pos <- acpData[cells,2,drop=TRUE];
      # Was: # pos <- acp[cells,2,drop=TRUE];

      # if samples to smooth, do smoothing
      if (nPos > 0) {
        matScorePos[[jj]] <- calcSmoothed(pos, dataList[[jj]][,sampPos,drop=FALSE], probeWindow=probeWindow, nProbes=nProbes, meanTrim=meanTrim);
      }
      if (nNeg > 0) {
        matScoreNeg[[jj]] <- calcSmoothed(pos, dataList[[jj]][,sampNeg,drop=FALSE], probeWindow=probeWindow, nProbes=nProbes, meanTrim=meanTrim);
      }
    } # for (jj ...)

    # Memory cleanup
    # Not needed anymore
    dataList <- NULL;

    verbose && exit(verbose);

    # MAT scores can be vectorized already here
    matScorePos <- unlist(matScorePos, use.names=FALSE);
    matScoreNeg <- unlist(matScoreNeg, use.names=FALSE);

    # calculate null distributions for scale factors
    if (length(sampNeg) > 0) {
      verbose && enter(verbose, "Calculating scale factor via null distributions");
      verbose && enter(verbose, "Gathering common info for calculating null distributions");
      chr <- acpData[allInds,1,drop=TRUE];
      pos <- acpData[allInds,2,drop=TRUE];
      # Faster than:
      ## chr <- acp[allInds,1, drop=TRUE];
      ## pos <- acp[allInds,2, drop=TRUE];
      verbose && exit(verbose);

      verbose && enter(verbose, "Calculating null distributions for controls and treatments in parallel");
      nullX <- cbind(neg=matScoreNeg, pos=matScorePos);
      nullDists <- calcNullDists(chr, pos, nullX);
      # Not needed anymore
      chr <- pos <- nullX <- NULL;
      verbose && str(verbose, nullDists);
      verbose && exit(verbose);

      scaleFactor <- nullDists$sd["pos"] / nullDists$sd["neg"];
      verbose && printf(verbose, "Scale factor: %.4g\n", scaleFactor);
      # Sanity check
      stopifnot(is.finite(scaleFactor));
      # Not needed anymore
      nullDists <- NULL;
      verbose && exit(verbose);
    } else {
      scaleFactor <- 1;
    } # if (length(sampNeg) > 0)


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Calculate MAT scores
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Calculating MAT scores (on intensity scale)");
    matScores <- matScorePos - scaleFactor*matScoreNeg;
    # Not needed anymore
    matScoreNeg <- matScorePos <- NULL;
    # ...on the intensity scale
    matScores <- 2^matScores;
    verbose && str(verbose, matScores);
    verbose && summary(verbose, matScores);
    verbose && exit(verbose);

    # Sanity check
    stopifnot(length(matScores) == length(allInds));

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Storing results
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Storing results");

    # Always validate the output file just before writing
    pathnameT <- Arguments$getWritablePathname(pathnameT, mustNotExist=!force);

    # Create CEL file to store results, if missing
    verbose && enter(verbose, "Creating CEL file for results, if missing");
    df <- getOneFile(ds);
    createFrom(df, filename=pathnameT, path=NULL, verbose=less(verbose));
    verbose && exit(verbose);

    verbose && enter(verbose, "Updating data file");
    verbose2 <- isVisible(verbose, -50);
    .updateCel(pathnameT, indices=allInds, intensities=matScores,
                                                        verbose=verbose2);
    verbose && exit(verbose);

    # Rename temporary file
    popTemporaryFile(pathnameT, verbose=verbose);

    ## Create checksum file
    dfZ <- getChecksumFile(pathname)

    # Not needed anymore
    filename <- pathnameT <- NULL;
    verbose && exit(verbose);


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Next column in design matrix
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Memory cleanup
    # Not needed anymore
    matScores <- NULL;

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);

    verbose && exit(verbose);
  } # for (ii ...)

  # Clean up
  # Not needed anymore
  cellsList <- allInds <- acpData <- NULL;

  outputDataSet <- getOutputDataSet(this, force=TRUE);

  # Sanity check
  stopifnot(length(outputDataSet) == ncol(design));

  verbose && exit(verbose);

  invisible(outputDataSet);
})


############################################################################
# HISTORY:
# 2012-09-02 [HB]
# o ROBUSTNESS: Added sanity check that process() for MatSmoothing returns
#   a dataset with the same number of arrays as number of columns in
#   the design matrix.
# 2009-06-28 [HB]
# o BUG FIX: Added missing getExpectedOutputFullnames() for MatSmoothing.
# 2009-05-27 [HB]
# o SPEED UP: Preloading acp[] data speeds things up > 5 times.
# o SPEED UP: Replaces two calcNullDist() calls with one calcNullDists()
#   call, which again is twice as fast.
# o SPEED UP: Made calcNullDist() twice as fast.
# o BUG FIX: calcNullDist() used the cell index instead of the cell
#   position for the first cell.
# 2009-05-27 [MR]
# o Added sanity check for CDF file that all units have exactly 1 group
# o Added 'stratifyBy="pm"' to the call to readUnits() and getCellIndices()
#   since the smoothing is only to be done on PM probes
# 2009-05-25 [HB]
# o Added getExpectedOutputFiles() special to MatSmoothing. Removed
#   isDone().
# 2009-05-23 [HB]
# o Added getOutputFiles() to MatSmoothing to work with new AromaTransform.
# o Now process() of MatSmoothing skips already process output files.
# o Updated process() of MatSmoothing to first write to a temporary file
#   which is then renamed.  This lower the risk for corrupt output files
#   due to processing interrupts.
# 2009-03-20 [MR]
# o Corrected the checking of isDone() for MatSmoothing objects.
# o If MatSmoothing has already been run, it returns the
#   AffymetrixCelSet object.
# 2009-01-13 [HB]
# o MEMORY CLEANUP: Cleaning out more "done" variables and earlier.
# o Code cleanup.
# 2008-03-21
# o Created from BackgroundCorrection.R.
############################################################################
