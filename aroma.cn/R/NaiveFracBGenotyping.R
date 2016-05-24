setConstructorS3("NaiveFracBGenotyping", function(dataSet=NULL, ...) {
  # Argument 'dataSet':
  if (!is.null(dataSet)) {
    dataSet <- Arguments$getInstanceOf(dataSet, "AromaUnitFracBCnBinarySet");
  }

  extend(Object(...), "NaiveFracBGenotyping",
    dataSet = dataSet
  )
})



setMethodS3("getDataSet", "NaiveFracBGenotyping", function(this, ...) {
  this$.ds;
})


setMethodS3("getGenotypeCallSet", "NaiveFracBGenotyping", function(this, ...) {
})


setMethodS3("nbrOfLoci", "NaiveFracBGenotyping", function(this, ...) {
  ds <- getDataSet(this);
  df <- getFile(ds, 1);
  nbrOfLoci(df);
})


setMethodS3("nbrOfSNPs", "NaiveFracBGenotyping", function(this, ...) {
  ds <- getDataSet(this);
  df <- getFile(ds, 1);
  nbrOfSNPs(df);
})


setMethodS3("processOne", "NaiveFracBGenotyping", function(this, df, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Calling genotypes");
  verbose && print(verbose, df);

  # To do
  gf <- NULL;

  verbose && print(verbose, gf);
  verbose && exit(verbose);

  gf;
})


setMethodS3("process", "NaiveFracBGenotyping", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Calling genotypes");
  ds <- getDataSet(this);
  verbose && cat(verbose, "Data set:");
  verbose && print(verbose, ds);

  tags <- "NGC";
  fullname <- paste(c(getFullName(ds), tags), collapse=",");
  chipType <- getChipType(ds, fullname=FALSE);
  outPath <- file.path("callData", fullname, chipType);
  outPath <- Arguments$getWritablePath(outPath);
  verbose && cat(verbose, "Output path: ", outPath);

  adjust <- 1.5;
  type <- NULL; rm(list="type"); # To please R CMD check

  platform <- getPlatform(df);
  chipType <- getChipType(df);
  nbrOfUnits <- nbrOfUnits(df);

  units <- NULL;
  for (kk in seq_along(ds)) {
    df <- getFile(ds, kk);
    verbose && enter(verbose, sprintf("Data file #%d ('%s') of %d", kk, getName(df), length(ds)));

    if (is.null(units)) {
      units <- seq_len(nbrOfUnits);

      # Identify units on ChrX and ChrY
      ugp <- getAromaUgpFile(ds);
      units23 <- getUnitsOnChromosome(ugp, 23);
      is23 <- is.element(units, units23);
      units24 <- getUnitsOnChromosome(ugp, 24);
      is24 <- is.element(units, units24);
    }

    # Create pathname of the resulting genotype call file
    tags <- getTags(df);
    tags <- setdiff(tags, "fracB");
    tags <- c(tags, "genotypes");
    fullname <- paste(c(getFullName(df), tags), collapse=",");
    filename <- sprintf("%s.acf", fullname);
    gcPathname <- Arguments$getWritablePathname(filename, path=outPath, mustNotExist=FALSE);

    # Create pathname of the resulting genotype confidence score file
    tags <- c(tags, "confidenceScores");
    fullname <- paste(c(getFullName(df), tags), collapse=",");
    filename <- sprintf("%s.acf", fullname);
    csPathname <- Arguments$getWritablePathname(filename, path=outPath, mustNotExist=FALSE);

    # Nothing to do?
    if (isFile(gcPathname) && isFile(csPathname)) {
      verbose && cat(verbose, "Already genotyped. Skipping.");
      verbose && exit(verbose);
      next;
    }

    verbose && enter(verbose, "Loading data");
    beta <- df[units,1,drop=TRUE];
    verbose && exit(verbose);

    # Call gender
    gender <- callXXorXY(beta[is23], beta[is24], adjust=adjust, from=0, to=1);

    # Call genotypes
    naValue <- as.double(NA);
    fit <- NULL;
    mu <- rep(naValue, times=length(units));
    cs <- rep(naValue, times=length(units));

    if (gender == "XY") {
      # All but ChrX & ChrY in male
      isDiploid <- (!(is23 | is24));
      use <- which(isDiploid);
      muT <- callNaiveGenotypes(beta[use], cn=2, adjust=adjust, from=0, to=1,
                                            verbose=less(verbose,10));
      fit <- attr(muT, 'modelFit');
      mu[use] <- muT;
      use <- which(!isDiploid);
      muT <- callNaiveGenotypes(beta[use], cn=1, adjust=adjust, from=0, to=1,
                                             verbose=less(verbose,10));
      mu[use] <- muT;
    } else {
      # All but ChrY in female
      isDiploid <- (!is24);
      use <- which(isDiploid);
      muT <- callNaiveGenotypes(beta[use], cn=2, adjust=adjust, from=0, to=1,
                                            verbose=less(verbose,10));
      fit <- attr(muT, 'modelFit');
      mu[use] <- muT;
    }

    verbose && print(verbose, table(mu, exclude=NULL));


    # Translate genotype calls in fracB space to (AA,AB,BB,...)
    calls <- rep(as.character(NA), times=length(mu));
    calls[mu ==   0] <- "AA";
    calls[mu == 1/2] <- "AB";
    calls[mu ==   1] <- "BB";
    verbose && print(verbose, table(calls, exclude=NULL));

    # Calculate confidence scores
    a <- fit$x[1];
    b <- fit$x[2];
    cs[isDiploid] <- rowMins(abs(cbind(beta[isDiploid]-a, beta[isDiploid]-b)));
    verbose && print(verbose, table(mu, exclude=NULL));


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Writing genotype calls (via temporary file)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Writing genotype calls");

    pathname <- gcPathname;
    pathnameT <- sprintf("%s.tmp", pathname);
    gf <- AromaUnitGenotypeCallFile$allocate(pathnameT, platform=platform, chipType=chipType, nbrOfRows=nbrOfUnits);
    footer <- readFooter(gf);
    footer$method <- "NaiveGenotypeCaller";
    writeFooter(gf, footer);
    # Not needed anymore
    footer <- NULL;

    updateGenotypes(gf, units=units, calls=calls);
    # Not needed anymore
    calls <- NULL;

    res <- file.rename(pathnameT, pathname);
    if (!isFile(pathname)) {
      throw("Failed to rename temporary file: ", pathnameT, " -> ", pathname);
    }
    if (isFile(pathnameT)) {
      throw("Failed to rename temporary file: ", pathnameT, " -> ", pathname);
    }
    # Not needed anymore
    pathnameT <- NULL

    gf <- AromaUnitGenotypeCallFile(pathname);
    print(verbose, gf);

    verbose && exit(verbose);


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Writing confidence scores (via temporary file)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Writing confidence scores");

    pathname <- csPathname;
    pathnameT <- sprintf("%s.tmp", pathname);
    csf <- AromaUnitSignalBinaryFile$allocate(pathnameT, platform=platform, chipType=chipType, nbrOfRows=nbrOfUnits, types="double", size=4, signed=TRUE);
    footer <- readFooter(csf);
    footer$method <- "NaiveGenotypeConfidenceScoreEstimator";
    writeFooter(csf, footer);
    # Not needed anymore
    footer <- NULL;

    csf[units, 1] <- cs;
    # Not needed anymore
    cs <- NULL;

    res <- file.rename(pathnameT, pathname);
    if (!isFile(pathname)) {
      throw("Failed to rename temporary file: ", pathnameT, " -> ", pathname);
    }
    if (isFile(pathnameT)) {
      throw("Failed to rename temporary file: ", pathnameT, " -> ", pathname);
    }
    # Not needed anymore
    pathnameT <- NULL;

    cf <- AromaUnitSignalBinaryFile(pathname);
    verbose && print(verbose, cf);
    verbose && exit(verbose);

    verbose && exit(verbose);
  }


  dataSet <- ""; # To do

  gcs <- AromaUnitGenotypeCallSet$byName(dataSet, tags="NGC", chipType="*");
  verbose && print(verbose, gcs);

  css <- AromaUnitSignalBinarySet$byName(dataSet, tags="NGC", chipType="*", pattern="confidenceScores", paths="callData");
  verbose && print(verbose, css);

  res <- list(calls=gcs, confidenceScores=css);

  verbose && exit(verbose);

  res;
})


############################################################################
# HISTORY:
# 2010-08-26
# o Added code to process().
# 2010-08-10
# o Created.
############################################################################
