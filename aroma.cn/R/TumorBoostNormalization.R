###########################################################################/**
# @RdocClass TumorBoostNormalization
#
# @title "The TumorBoostNormalization class"
#
# \description{
#  @classhierarchy
#
#  TumorBoost is normalization method that normalizes the allele B fractions
#  of a tumor sample given the allele B fractions and genotype calls for
#  a matched normal.
#  The method is a single-sample (single-pair) method.  It does not require
#  total copy number estimates.
#  The normalization is done such that the total copy number is unchanged
#  afterwards.
# }
#
# @synopsis
#
# \arguments{
#  \item{dsT}{An @see "aroma.core::AromaUnitFracBCnBinarySet" of
#     tumor samples.}
#  \item{dsN}{An @see "aroma.core::AromaUnitFracBCnBinarySet" of
#     match normal samples.}
#  \item{gcN}{An @see "aroma.core::AromaUnitGenotypeCallSet" of
#     genotypes for the normals.}
#  \item{flavor}{A @character string specifying the type of
#     correction applied.}
#  \item{preserveScale}{If @TRUE, SNPs that are heterozygous in the
#    matched normal are corrected for signal compression using an estimate
#    of signal compression based on the amount of correction performed
#    by TumorBoost on SNPs that are homozygous in the matched normal.}
#  \item{collapseHomozygous}{If @TRUE, SNPs that are homozygous in the
#    matched normal are also called homozygous in the tumor, that is,
#    it's allele B fraction is collapsed to either 0 or 1.
#    If @FALSE, the homozygous values are normalized according the
#    model. [NOT USED YET]
#  }
#  \item{tags}{(Optional) Sets the tags for the output data sets.}
#  \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB, PN"
#*/###########################################################################
setConstructorS3("TumorBoostNormalization", function(dsT=NULL, dsN=NULL, gcN=NULL, flavor=c("v4", "v3", "v2", "v1"), preserveScale=TRUE, collapseHomozygous=FALSE, tags="*", ...) {
  # Validate arguments
  if (!is.null(dsT)) {
    # Argument 'flavor':
    flavor <- match.arg(flavor);

    # Arguments 'dsT' and 'dsN'
    dsList <- list(dsT=dsT, dsN=dsN);
    className <- "AromaUnitFracBCnBinarySet";
    for (kk in seq_along(dsList)) {
      key <- names(dsList)[kk];
      ds <- dsList[[kk]];
      ds <- Arguments$getInstanceOf(ds, className, .name="key");
    }

    # Argument 'gcN':
    gcN <- Arguments$getInstanceOf(gcN, "AromaUnitGenotypeCallSet");

    # Assert that each data set contains the same number of files
    dsList$gcN <- gcN;
    for (jj in 1:(length(dsList)-1L)) {
      keyJJ <- names(dsList)[jj];
      dsJJ <- dsList[[jj]];
      nJJ <- length(dsJJ);
      chipTypeJJ <- getChipType(dsJJ);
      for (kk in (jj+1):length(dsList)) {
        keyKK <- names(dsList)[kk];
        dsKK <- dsList[[kk]];
        nKK <- length(dsKK);
        chipTypeKK <- getChipType(dsKK);

        # Assert that each data set contains the same number of files
        if (nKK != nJJ) {
          throw(sprintf("The number of files in '%s' and '%s' does not match: %s != %s", keyKK, keyJJ, nKK, nJJ));
        }

        # Assert that each data set is for the same chip type
        if (chipTypeKK != chipTypeJJ) {
          throw(sprintf("The chip types for '%s' and '%s' does not match: %s != %s", keyKK, keyJJ, chipTypeKK, chipTypeJJ));
        }
      } # for (kk ...)
    } # for (jj ...)
  } # if (!is.null(dsT))

  preserveScale <- Arguments$getLogical(preserveScale);

  collapseHomozygous <- Arguments$getLogical(collapseHomozygous);
  if (collapseHomozygous) {
    throw("collapseHomozygous=FALSE is currently not implemented.");
  }

  # Arguments '...':
  args <- list(...);
  if (length(args) > 0L) {
    argsStr <- paste(names(args), collapse=", ");
    throw("Unknown arguments: ", argsStr);
  }

  this <- extend(Object(...), "TumorBoostNormalization",
    .dsT = dsT,
    .dsN = dsN,
    .gcN = gcN,
    .flavor = flavor,
    .preserveScale = preserveScale
  );

  setTags(this, tags);

  this;
})


setMethodS3("as.character", "TumorBoostNormalization", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1L]);

  s <- c(s, sprintf("Flavor: %s", getFlavor(this)));

  dsList <- getDataSets(this);
  s <- c(s, sprintf("Data sets (%d):", length(dsList)));
  for (kk in seq_along(dsList)) {
    ds <- dsList[[kk]];
    s <- c(s, sprintf("<%s>:", capitalize(names(dsList)[kk])));
    s <- c(s, as.character(ds));
  }

  GenericSummary(s);
}, protected=TRUE)


setMethodS3("getAsteriskTags", "TumorBoostNormalization", function(this, collapse=NULL, ...) {
  tags <- "TBN";

  flavor <- getFlavor(this);
  if (flavor != "v4") {
    tags <- c(tags, flavor);
  }

  preserveScale <- this$.preserveScale;
  if (!preserveScale) {
    tags <- c(tags, "ns");
  }

  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  }

  tags;
}, protected=TRUE)


setMethodS3("getName", "TumorBoostNormalization", function(this, ...) {
  ds <- getInputDataSet(this);
  getName(ds);
})

setMethodS3("getFlavor", "TumorBoostNormalization", function(this, ...) {
  this$.flavor;
}, protected=TRUE)


setMethodS3("getTags", "TumorBoostNormalization", function(this, collapse=NULL, ...) {
  # "Pass down" tags from input data set
  ds <- getInputDataSet(this);
  tags <- getTags(ds, collapse=collapse);

  # Get class-specific tags
  tags <- c(tags, this$.tags);

  # Update default tags
  tags[tags == "*"] <- getAsteriskTags(this, collapse=",");

  # Collapsed or split?
  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  } else {
    tags <- unlist(strsplit(tags, split=","));
  }

  if (length(tags) == 0L)
    tags <- NULL;

  tags;
})


setMethodS3("setTags", "TumorBoostNormalization", function(this, tags="*", ...) {
  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));
    tags <- tags[nchar(tags) > 0L];
  }

  this$.tags <- tags;
})


setMethodS3("getFullName", "TumorBoostNormalization", function(this, ...) {
  name <- getName(this);
  tags <- getTags(this);
  fullname <- paste(c(name, tags), collapse=",");
  fullname <- gsub("[,]$", "", fullname);
  fullname;
})


setMethodS3("getDataSets", "TumorBoostNormalization", function(this, ...) {
  list(tumor=this$.dsT, normal=this$.dsN, normalCalls=this$.gcN);
}, protected=TRUE)

setMethodS3("getInputDataSet", "TumorBoostNormalization", function(this, ...) {
  this$.dsT;
})

setMethodS3("getNormalDataSet", "TumorBoostNormalization", function(this, ...) {
  this$.dsN;
})

setMethodS3("getNormalGenotypeCallSet", "TumorBoostNormalization", function(this, ...) {
  this$.gcN;
})

setMethodS3("getRootPath", "TumorBoostNormalization", function(this, ...) {
  "totalAndFracBData";
}, protected=TRUE)


setMethodS3("getPath", "TumorBoostNormalization", function(this, create=TRUE, ...) {
  # Create the (sub-)directory tree for the data set

  # Root path
  rootPath <- getRootPath(this);

  # Full name
  fullname <- getFullName(this);

  # Chip type
  ds <- getInputDataSet(this);
  chipType <- getChipType(ds, fullname=FALSE);

  # The full path
  path <- filePath(rootPath, fullname, chipType);

  # Create path?
  if (create) {
    path <- Arguments$getWritablePath(path);
  } else {
    path <- Arguments$getReadablePath(path, mustExist=FALSE);
  }

  # Verify that it is not the same as the input path
  inPath <- getPath(getInputDataSet(this));
  if (getAbsolutePath(path) == getAbsolutePath(inPath)) {
    throw("The generated output data path equals the input data path: ", path, " == ", inPath);
  }

  path;
}, protected=TRUE)


setMethodS3("nbrOfFiles", "TumorBoostNormalization", function(this, ...) {
  ds <- getInputDataSet(this);
  length(ds);
})


setMethodS3("getOutputDataSet", "TumorBoostNormalization", function(this, ...) {
  ds <- getInputDataSet(this);
  path <- getPath(this);
  res <- byPath(ds, path=path, ...);
  res;
})


setMethodS3("process", "TumorBoostNormalization", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  units <- NULL;

  verbose && enter(verbose, "TumorBoost normalization");
  nbrOfFiles <- nbrOfFiles(this);
  verbose && cat(verbose, "Number of arrays: ", nbrOfFiles);

  flavor <- getFlavor(this);
  verbose && cat(verbose, "Flavor: ", flavor);

  dsList <- getDataSets(this);
  chipType <- getChipType(dsList[[1L]], fullname=FALSE);
  verbose && cat(verbose, "Chip type: ", chipType);

  outPath <- getPath(this);
  for (kk in seq_len(nbrOfFiles)) {
    dfList <- lapply(dsList, FUN=getFile, kk);
    dfT <- dfList$tumor;
    name <- getName(dfT);
    verbose && enter(verbose, sprintf("Sample #%d ('%s') of %d",
                                                    kk, name, nbrOfFiles));

    # Output file
    filename <- getFilename(dfT);
    pathname <- Arguments$getReadablePathname(filename, path=outPath, mustExist=FALSE);

    # Nothing to do?
    if (isFile(pathname)) {
      verbose && cat(verbose, "Already normalized.");
      verbose && exit(verbose);
      next;
    }

    verbose && print(verbose, dfList);

    if (is.null(units)) {
      verbose && enter(verbose, "Identifying units to read");
      units <- seq_len(nbrOfUnits(dfList$tumor));
      verbose && exit(verbose);
    }

    verbose && enter(verbose, "Reading all data");
    betaT <- dfList$tumor[units,1L,drop=TRUE];
    keep <- is.finite(betaT);
    unitsT <- units[keep];
    betaT <- betaT[keep];
    verbose && cat(verbose, "Allele B fractions for the tumor:");
    verbose && str(verbose, betaT);

    verbose && cat(verbose, "Allele B fractions for the matched normal:");
    betaN <- dfList$normal[unitsT,1,drop=TRUE];
    verbose && str(verbose, betaN);

    verbose && cat(verbose, "Genotypes for the matched normal:");
    gfN <- dfList$normalCalls;
    muN <- extractGenotypes(gfN, units=unitsT, encoding="fracB", drop=TRUE);
    verbose && str(verbose, muN);
    verbose && print(verbose, table(muN));
    verbose && exit(verbose);


    verbose && enter(verbose, "Normalizing tumor allele B fractions");
    verbose && cat(verbose, "Flavor: ", flavor);
    verbose && enter(verbose, "Estimating SNP effects");
    # NOTE: It is possible that 'delta' has NA:s.
    delta <- (betaN - muN);
    verbose && str(verbose, delta);
    verbose && exit(verbose);

    verbose && enter(verbose, "Rescaling correction factor");
    if (flavor == "v1") {
      b <- 1;
    } else if (flavor == "v2") {
      b <- rep(1, times=length(delta));
      isDown <- (betaT < betaN);
      idxs <- which(isDown);
      # NOTE: It is possible that 'b' has NA:s.
      b[idxs] <- betaT[idxs]/betaN[idxs];
      idxs <- which(!isDown);
      # NOTE: It is possible that 'b' has NA:s.
      b[idxs] <- (1-betaT[idxs])/(1-betaN[idxs]);
      # Not needed anymore
      isDown <- isHomA <- isHomB <- idxs <- NULL;
    } else if (flavor == "v3") {
      b <- rep(1, times=length(delta));
      isHomA <- (muN == 0);
      isHomB <- (muN == 1);
      isHet <- !isHomA & !isHomB;
      isDown <- (betaT < betaN);
      idxs <- which((isHet & isDown) | isHomA);
      # NOTE: It is possible that 'b' has NA:s.
      b[idxs] <- betaT[idxs]/betaN[idxs];
      idxs <- which((isHet & !isDown) | isHomB);
      # NOTE: It is possible that 'b' has NA:s.
      b[idxs] <- (1-betaT[idxs])/(1-betaN[idxs]);
      # Not needed anymore
      isDown <- isHet <- isHomA <- isHomB <- idxs <- NULL;
    } else if (flavor == "v4") {
      b <- rep(1, times=length(delta));
      isHet <- (muN != 0 & muN != 1);
      isDown <- (betaT < betaN);
      idxs <- which(isHet & isDown);
      # NOTE: It is possible that 'b' has NA:s.
      b[idxs] <- betaT[idxs]/betaN[idxs];
      idxs <- which(isHet & !isDown);
      # NOTE: It is possible that 'b' has NA:s.
      b[idxs] <- (1-betaT[idxs])/(1-betaN[idxs]);
      # Not needed anymore
      isDown <- isHet <- idxs <- NULL;
    }
    verbose && cat(verbose, "Scaling factor:");
    verbose && str(verbose, b);
    verbose && summary(verbose, b);
    verbose && exit(verbose);

    verbose && enter(verbose, "Normalizing");
    # NOTE: It is possible that we introduce NA:s via 'b' and 'delta'.
    betaTN <- betaT - b*delta;

    preserveScale <- this$.preserveScale;
    if (preserveScale) {
      verbose && enter(verbose, "Correcting for signal compression");

      isHom <- (muN == 0 | muN == 1);
      idxs <- which(isHom);
      eta <- median(abs(betaT[idxs]-1/2), na.rm=TRUE);
      verbose && cat(verbose, "Signal compression in homozygous SNPs before TBN");
      verbose && str(verbose, 1/2-eta);
      etaC <- median(abs(betaTN[idxs]-1/2), na.rm=TRUE);
      verbose && cat(verbose, "Signal compression in homozygous SNPs after TBN");
      verbose && str(verbose, 1/2-etaC);

      # Correction factor
      sf <- etaC/eta;

      isHet <- !isHom;
      isDown <- (betaTN < 1/2);
      idxs <- which(isHet & isDown);
      betaTN[idxs] <- 1/2 - sf * (1/2 - betaTN[idxs]);
      idxs <- which(isHet & !isDown);
      betaTN[idxs] <- 1/2 + sf * (betaTN[idxs] - 1/2);

      # Not needed anymore
      isDown <- isHom <- isHet <- idxs <- eta <- etaC <- sf <- NULL;
      verbose && exit(verbose);
    }
    verbose && str(verbose, betaTN);
    verbose && exit(verbose);
    verbose && exit(verbose);


    verbose && enter(verbose, "Storing normalized data");

    verbose && enter(verbose, "Allocating to temporary file");
    pathnameT <- sprintf("%s.tmp", pathname);
    dfT <- dfList$tumor;
    outPath <- Arguments$getWritablePath(outPath);
    file.copy(getPathname(dfT), pathnameT);
    dfTC <- newInstance(dfT, pathnameT);
    srcFiles <- lapply(dfList, FUN=function(df) {
      list(
        filename = getFilename(df),
        filesize = getFileSize(df),
        checksum = getChecksum(df)
      )
    });
    footer <- readFooter(dfTC);
    footer$srcFiles <- footer;
    writeFooter(dfTC, footer);
    verbose && exit(verbose);

    verbose && enter(verbose, "Writing to temporary file");
    dfTC[unitsT,1] <- betaTN;
    verbose && exit(verbose);

    # Renaming
    verbose && enter(verbose, "Renaming temporary file");
    file.rename(pathnameT, pathname);
    if (!isFile(pathname)) {
      throw("Failed to rename temporary file: ", pathnameT, " -> ", pathname);
    }
    verbose && exit(verbose);
    verbose && exit(verbose);

#    verbose && enter(verbose, "Storing estimates to priorData/");
#    path <- file.path("priorData", "chipTypes", chipType);
#    path <- Arguments$getWritablePath(path);
#    verbose && exit(verbose);

    verbose && exit(verbose);
  } # for (kk ...)

  res <- getOutputDataSet(this, verbose=less(verbose, 1));

  verbose && exit(verbose);

  invisible(res);
})

############################################################################
# HISTORY:
# 2010-08-04 [PN]
# o ROBUSTNESS: Added 'na.rm=TRUE' to 'median'.
# o CLEAN UP: Removed an unnecessary 'rm'.
# o Added option 'preserveScale' to correct for signal compression in
#   heterozygous SNPs.  Defaults to 'TRUE'.
# 2010-06-20
# o CLEAN UP: Removed a duplicated line of code.
# 2009-12-09
# o Made flavor="v4" of TumorBoostNormalization the default, and if used
#   then no "flavor" tag is added.
# 2009-09-11
# o Added table verbose output of the read genotypes.  This is just so one
#   can verify the encoding, i.e. 0, 1/2, or 1.
# 2009-07-15
# o BUG FIX: TumorBoostNormalization: the 'srcFiles' attribute in file
#   footer of the result files contained a duplicated default footer
#   instead of the tumor-normal pair.
# 2009-07-02
# o Added model 'flavor' "v4" which corrects heterozygots according to "v2"
#   and homozygotes according to "v1".
# o Added model 'flavor' "v3".  Suggested by PN last night over a Guinness
#   at the pub after a long day of hard work.
# 2009-06-22
# o Added model 'flavor' "v2".
# 2009-06-08
# o The constructor of TumorBoostNormalization now only takes an
#   AromaUnitGenotypeCallSet for argument 'gcN'.  It no longer takes an
#   AromaUnitFracBCnBinarySet object.
# 2009-05-17
# o Now the constructor of TumorBoostNormalization asserts that there are
#   no stray arguments.
# 2009-04-29
# o Created.
############################################################################
