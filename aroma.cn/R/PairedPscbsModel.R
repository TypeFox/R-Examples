###########################################################################/**
# @RdocClass PairedPscbsModel
#
# @title "The PairedPscbsModel class"
#
# \description{
#  @classhierarchy
#
#  This class represents the Paired PSCBS method [1], which
#  segments matched tumor-normal parental copy-number data into
#  piecewise constant segments.
# }
#
# @synopsis
#
# \arguments{
#   \item{dsT, dsN}{The tumor and the normal
#     @see "aroma.core::AromaUnitPscnBinarySet".}
#   \item{tags}{Tags added to the output data sets.}
#   \item{...}{(Optional) Additional arguments passed to
#     @see "PSCBS::segmentByPairedPSCBS".}
#   \item{dropTcnOutliers}{If @TRUE, then TCN outliers are dropped using
#     @see "PSCBS::dropSegmentationOutliers".}
#   \item{gapMinLength}{Genomic regions with no data points that are of
#     this length and greater are considered to be "gaps" and are ignored
#     in the segmentation.  If +@Inf, no gaps are identified.}
#   \item{seed}{An optional @integer specifying the random seed to be
#     used in the segmentation.  Seed needs to be set for exact numerical
#     reproducibility.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \examples{\dontrun{
#   @include "../incl/PairedPscbsModel.Rex"
# }}
#
# \references{
#  [1] ... \cr
# }
#
# \seealso{
#   ...
# }
#
#*/###########################################################################
setConstructorS3("PairedPscbsModel", function(dsT=NULL, dsN=NULL, tags="*", ..., dropTcnOutliers=TRUE, gapMinLength=1e6, seed=NULL) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Load required packages
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(dsT)) {
    .requirePkg("PSCBS", quietly=TRUE);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(dsT)) {
    # Argument 'dsT':
    dsT <- Arguments$getInstanceOf(dsT, "AromaUnitPscnBinarySet");

    # Argument 'dsN':
    dsN <- Arguments$getInstanceOf(dsN, "AromaUnitPscnBinarySet");

    # Assert that the chip types are compatile
    if (getChipType(dsT) != getChipType(dsN)) {
        throw("Argument 'dsT' and 'dsN' are of different chip types: ",
              getChipType(dsT), " != ", getChipType(dsN));
    }

    # Assert that the ds sets have the same number ds files
    nbrOfFiles <- length(dsT);
    if (nbrOfFiles != length(dsN)) {
      throw("The number of samples in 'dsT' and 'dsN' differ: ",
             nbrOfFiles, " != ", length(dsN));
    }
  }

  # Argument 'dropTcnOutliers':
  dropTcnOutliers <- Arguments$getLogical(dropTcnOutliers);

  # Argument 'gapMinLength':
  gapMinLength <- Arguments$getDouble(gapMinLength, range=c(0,Inf));

  # Argument 'seed':
  if (!is.null(seed)) {
    seed <- Arguments$getInteger(seed);
  }

  # Arguments '...'; optional arguments to segmentByPairedPSCBS()
  extraArgs <- list(...);
  if (length(extraArgs) > 0) {
    keys <- names(extraArgs);
    if (is.null(keys)) {
      throw("Optional arguments to PairedPscbsModel passed via '...' must be named.");
    }

    nok <- which(nchar(keys) == 0);
    if (length(nok) > 0) {
      throw("All arguments to PairedPscbsModel passed via '...' must be named.");
    }
  }

  this <- extend(Object(), c("PairedPscbsModel", uses("ParametersInterface")),
    .dsT = dsT,
    .dsN = dsN,
    .extraArgs = extraArgs,
    .dropTcnOutliers = dropTcnOutliers,
    .gapMinLength = gapMinLength,
    .seed = seed
  );

  setTags(this, tags);

  this;
})


setMethodS3("as.character", "PairedPscbsModel", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);

  dsList <- getDataSets(this);
  for (key in names(dsList)) {
    ds <- dsList[[key]];
    s <- c(s, sprintf("%s data set:", capitalize(key)));
    s <- c(s, as.character(ds));
  }

  nbrOfFiles <- nbrOfFiles(this);
  s <- c(s, sprintf("Number of arrays: %d", nbrOfFiles));
   s <- c(s, sprintf("Additional parameters: %s", getParametersAsString(this)));

  GenericSummary(s);
}, protected=TRUE)



setMethodS3("getRandomSeed", "PairedPscbsModel", function(this, ...) {
  this$.seed;
}, protected=TRUE)

setMethodS3("setRandomSeed", "PairedPscbsModel", function(this, seed, ...) {
  # Argument 'seed':
  if (!is.null(seed)) {
    seed <- Arguments$getInteger(seed);
  }

  this$.seed <- seed;
  invisible(this);
}, protected=TRUE)


setMethodS3("getTumorDataSet", "PairedPscbsModel", function(this, ...) {
  this$.dsT;
})

setMethodS3("getNormalDataSet", "PairedPscbsModel", function(this, ...) {
  this$.dsN;
})

setMethodS3("getDataSets", "PairedPscbsModel", function(this, ...) {
  list(tumor=getTumorDataSet(this), normal=getNormalDataSet(this));
})


setMethodS3("getOutputDataSet", "PairedPscbsModel", function(this, ...) {
  path <- getPath(this);
  res <- PairedPSCBSFileSet$byPath(path);
  res;
})


setMethodS3("getFitFunction", "PairedPscbsModel", function(this, ...) {
  use("PSCBS (>= 0.43.0)")
  segmentByPairedPSCBS <- PSCBS::segmentByPairedPSCBS

  defaultSeed <- getRandomSeed(this);
  fitFcn <- function(..., seed=defaultSeed) {
    segmentByPairedPSCBS(..., seed=seed);
  }
  fitFcn;
}, protected=TRUE)

setMethodS3("getAsteriskTags", "PairedPscbsModel", function(this, collapse=NULL, ...) {
  tags <- "PairedPSCBS";
  params <- getParameters(this);
  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  }

  tags;
}, protected=TRUE)


setMethodS3("getName", "PairedPscbsModel", function(this, ...) {
  dsT <- getTumorDataSet(this);
  getName(dsT);
})



setMethodS3("getTags", "PairedPscbsModel", function(this, collapse=NULL, ...) {
  # "Pass down" tags from input tumor data set
  dsT <- getTumorDataSet(this);
  tags <- getTags(dsT, collapse=collapse);

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

  if (length(tags) == 0)
    tags <- NULL;

  tags;
})


setMethodS3("setTags", "PairedPscbsModel", function(this, tags="*", ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));
    tags <- tags[nchar(tags) > 0];
  }

  this$.tags <- tags;
})


setMethodS3("getFullName", "PairedPscbsModel", function(this, ...) {
  name <- getName(this);
  tags <- getTags(this);
  fullname <- paste(c(name, tags), collapse=",");
  fullname <- gsub("[,]$", "", fullname);
  fullname;
})

setMethodS3("getChipType", "PairedPscbsModel", function(this, ...) {
  dsT <- getTumorDataSet(this);
  getChipType(dsT, ...);
})


setMethodS3("nbrOfFiles", "PairedPscbsModel", function(this, ...) {
  dsT <- getTumorDataSet(this);
  length(dsT);
})



setMethodS3("getRootPath", "PairedPscbsModel", function(this, ...) {
  "pscbsData";
}, protected=TRUE)


setMethodS3("getPath", "PairedPscbsModel", function(this, create=TRUE, ...) {
  # Create the (sub-)directory tree for the data set

  # Root path
  rootPath <- getRootPath(this);

  # Full name
  fullname <- getFullName(this);

  # Chip type
  chipType <- getChipType(this, fullname=FALSE);

  # The full path
  path <- filePath(rootPath, fullname, chipType);

  # Create path?
  if (create) {
    path <- Arguments$getWritablePath(path);
  } else {
    path <- Arguments$getReadablePath(path, mustExist=FALSE);
  }

  # Verify that it is not the same as the input path
  dsList <- getDataSets(this);
  inPath <- getPath(dsList$tumor);
  if (getAbsolutePath(path) == getAbsolutePath(inPath)) {
    throw("The generated output data path equals the input data path: ", path, " == ", inPath);
  }

  path;
}, protected=TRUE)



setMethodS3("getParameters", "PairedPscbsModel", function(this, ...) {
  params <- NextMethod("getParameters");
  params$dropTcnOutliers <- this$.dropTcnOutliers;
  params$gapMinLength <- this$.gapMinLength;
  params$seed <- getRandomSeed(this);
  params <- c(params, this$.extraArgs);
  params;
}, protected=TRUE);


setMethodS3("getOptionalArguments", "PairedPscbsModel", function(this, ...) {
  args <- getParameters(this);
  args$dropTcnOutliers <- NULL;
  args$gapMinLength <- NULL;
  args;
}, protected=TRUE);



setMethodS3("getChromosomes", "PairedPscbsModel", function(this, ...) {
  dsT <- getTumorDataSet(this);
  ugp <- getAromaUgpFile(dsT);
  getChromosomes(ugp, ...);
})


setMethodS3("indexOf", "PairedPscbsModel", function(this, ...) {
  dsT <- getTumorDataSet(this);
  indexOf(dsT, ...);
})


setMethodS3("fit", "PairedPscbsModel", function(this, arrays=NULL, chromosomes=getChromosomes(this), maxNAFraction=1/5, force=FALSE, ..., .retResults=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'arrays':
  if (identical(arrays, "fitted")) {
  } else {
    arrays <- indexOf(this, arrays);
  }

  allChromosomes <- getChromosomes(this);

  # Argument 'chromosomes':
  if (identical(chromosomes, "fitted")) {
  } else if (is.null(chromosomes)) {
    chromosomes <- getChromosomes(this);
  } else if (is.numeric(chromosomes)) {
    chromosomes <- Arguments$getChromosomes(chromosomes,
                                                range=range(allChromosomes));
##    chromosomes <- as.character(chromosomes);  ## TODO
##    chromosomes[chromosomes == "23"] <- "X";   ## TODO
    chromosomes <- intersect(chromosomes, allChromosomes);
  } else if (is.character(chromosomes)) {
    chromosomes <- Arguments$getChromosomes(chromosomes,
                                                range=range(allChromosomes));
##    chromosomes[chromosomes == "23"] <- "X";   ## TODO
    chromosomes <- intersect(chromosomes, getChromosomes(this));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  path <- getPath(this);
  mkdirs(path);

  fitFcn <- getFitFunction(this, verbose=less(verbose, 50));

  params <- getParameters(this);
  verbose && cat(verbose, "Parameters:");
  verbose && str(verbose, params);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieving tumor-normal data sets
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dsList <- getDataSets(this, verbose=verbose);
  verbose && cat(verbose, "Data sets:");
  verbose && print(verbose, dsList);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract annotation data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ugp <- getAromaUgpFile(dsList$tumor);
  cp <- readDataFrame(ugp, verbose=less(verbose, 50));
#  colnames(cp) <- gsub("position", "x", colnames(cp), fixed=TRUE);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Array by array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res <- list();
  nbrOfArrays <- length(arrays);
  for (aa in seq_len(nbrOfArrays)) {
    array <- arrays[aa];

    # The tumor and normal data files
    dfT <- getFile(dsList$tumor, array);
    dfN <- getFile(dsList$normal, array);

    names <- c(T=getName(dfT), N=getName(dfN));
    pairName <- paste(unique(names), collapse="_vs_");

    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d",
                                                aa, pairName, nbrOfArrays));

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Generate fullname
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Paired tags
    # Generate pair tags (e.g. "abc_vs_def")
    tagsT <- setdiff(getTags(dfT), c("pscn", "T", "N"));
    tagsN <- setdiff(getTags(dfN), c("pscn", "T", "N"));
    tagsT <- paste(tagsT, collapse="-");
    tagsN <- paste(tagsN, collapse="-");
    pairTags <- paste(c(tagsT, tagsN), collapse="_vs_");
    pairTags <- c("TvsN", pairTags);
    fullname <- paste(c(pairName, pairTags), collapse=",");

    verbose && cat(verbose, "Tags: ", pairTags);
    verbose && cat(verbose, "Fullname: ", fullname);


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Get pathname
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    filename <- sprintf("%s,PairedPSCBS.xdr", fullname);
    pathname <- filePath(path, filename);

    # Already done?
    if (!force && isFile(pathname)) {
      verbose && enter(verbose, "Loading results from file");
      verbose && cat(verbose, "Pathname: ", pathname);
      fit <- loadObject(pathname);
      verbose && cat(verbose, "Fit object: ", class(fit)[1]);
      verbose && exit(verbose);
    } else {
      # Time the fitting.
      startTime <- processTime();
      timers <- list(total=0, read=0, fit=0, write=0, gc=0);
      tTotal <- processTime();

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Extract the tumor-normal data
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Extracting tumor-normal data");

      tRead <- processTime();
      dataT <- readDataFrame(dfT, verbose=less(verbose,50));
      dataN <- readDataFrame(dfN, verbose=less(verbose,50));
      timers$read <- timers$read + (processTime() - tRead);

      data <- data.frame(chromosome=cp$chromosome, x=cp$position,
                         CT=2*dataT$total/dataN$total,
                         betaT=dataT$fracB, betaN=dataN$fracB);
      verbose && str(verbose, data);

      # Not needed anymore
      dfT <- dfN <- dataT <- dataN <- NULL;
      verbose && exit(verbose);


      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Drop TCN outliers
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      dropTcnOutliers <- this$.dropTcnOutliers;
      if (dropTcnOutliers) {
        verbose && enter(verbose, "Dropping TCN outliers");
        data <- dropSegmentationOutliers(data);
        verbose && exit(verbose);
      }

      nbrOfLoci <- nrow(data);
      verbose && cat(verbose, "Number of loci: ", nbrOfLoci);


      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Find large gaps
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      minLength <- params$gapMinLength;
      if (is.finite(minLength)) {
        verbose && enter(verbose, "Finding large gaps");
        verbose && cat(verbose, "Minimal gap length: ", minLength);
        gaps <- findLargeGaps(data, minLength=minLength);
        verbose && print(verbose, gaps);
        knownSegments <- gapsToSegments(gaps, dropGaps=FALSE);
        # Not needed anymore
        gaps <- NULL;
        verbose && exit(verbose);
      }

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Fit segmentation model
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Calling model fit function");

      optArgs <- getOptionalArguments(this);
      verbose && cat(verbose, "Optional arguments (may be ignored/may give an error/warning):");
      verbose && str(verbose, optArgs);
      args <- list(data, knownSegments=knownSegments);
      args <- c(args, optArgs);
      args <- c(args, list(...));
      verbose && cat(verbose, "All arguments:");
      verbose && str(verbose, args);
      args <- c(args, list(...), list(verbose=less(verbose, 1)));
      tFit <- processTime();
      fit <- do.call("fitFcn", args);
      fit <- setSampleName(fit, pairName);
      verbose && str(verbose, fit);
      timers$fit <- timers$fit + (processTime() - tFit);
      # Not needed anymore
      data <- NULL;

      verbose && cat(verbose, "Class of fitted object: ", class(fit)[1]);
      verbose && printf(verbose, "Time to fit segmentation model: %.2fmin\n", timers$fit[3]/60);

      verbose && exit(verbose);


      # Garbage collection
      tGc <- processTime();
      gc <- gc();
      timers$gc <- timers$gc + (processTime() - tGc);
      verbose && print(verbose, gc);

      verbose && enter(verbose, "Saving to file");
      verbose && cat(verbose, "Pathname: ", pathname);
      tWrite <- processTime();
      saveObject(fit, file=pathname);
      timers$write <- timers$write + (processTime() - tWrite);
      verbose && exit(verbose);

      timers$total <- timers$total + (processTime() - tTotal);

      # Report time profiling
      totalTime <- processTime() - startTime;

      if (verbose) {
        t <- totalTime[3];
        printf(verbose, "Total time: %.2fs == %.2fmin\n", t, t/60);
        t <- totalTime[3]/nbrOfLoci;
        printf(verbose, "Total time per 1000 locus (with %d loci): %.2fs\n", nbrOfLoci, 1000*t);
        # Get distribution of what is spend where
        t <- lapply(timers, FUN=function(timer) unname(timer[3]));
        t <- unlist(t);
        t <- 100 * t / t["total"];
        printf(verbose, "Fraction of time spent on different tasks: Fitting: %.1f%%, Reading: %.1f%%, Writing: %.1f%%, Explicit garbage collection: %.1f%%\n", t["fit"], t["read"], t["write"], t["gc"]);
      }
    } # ...

    hookName <- "onFit.PairedPscbsModel";
    verbose && enter(verbose, sprintf("Calling %s() hooks", hookName));
    callHooks(hookName, fit=fit, fullname=fullname);
    verbose && exit(verbose);

    if (.retResults) {
      res[[pairName]] <- fit;
      # Not needed anymore
      fit <- NULL;
    }

    verbose && exit(verbose);
  } # for (aa in ...)

  invisible(res);
}) # fit()


############################################################################
# HISTORY:
# 2013-08-12
# o BUG FIX: getPath() for PairedPscbsModel would throw an error on
#   getInputDataSet() not defined.
# 2012-11-21
# o Now class utilizes the new ParametersInterface.
# 2012-09-20
# o Now PairedPscbsModel() excludes the actual gaps from the known
#   segments it passes to segmentByPairedPSCBS().
# 2012-09-19
# o Now getOutputDataSet() for PairedPscbsModel returns a
#   PairedPSCBSFileSet.
# 2012-09-15
# o Now fit() for PairedPscbsModel generates pair names iff tumor and
#   normal names don't match, e.g. 'GSM517071_vs_GSM517072' (if match
#   then just 'Patient1').  It also generated "pair" tags.
# 2012-07-22
# o fit() seems to work now.
# o Now utilizing the new AromaUnitPscnBinarySet class.
# 2012-07-20
# o Created from CalMaTeModel.
############################################################################
