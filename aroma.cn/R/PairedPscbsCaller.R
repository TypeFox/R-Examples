setConstructorS3("PairedPscbsCaller", function(dataSet=NULL, tags="*", calls=c("ROH", "AB", "LOH"), ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Load required packages
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(dataSet)) {
    .requirePkg("PSCBS", quietly=TRUE);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  if (!is.null(dataSet)) {
    dataSet <- Arguments$getInstanceOf(dataSet, "PairedPSCBSFileSet");
  }

  # Argument 'calls':
  calls <- match.arg(calls, several.ok=TRUE);

  # Arguments '...':
  optionalArgs <- list();

  extend(AromaTransform(dataSet=dataSet, tags=tags,
               .reqSetClass="PairedPSCBSFileSet"), "PairedPscbsCaller",
    .calls = calls,
    .optionalArgs = optionalArgs
  );
}) # PairedPscbsCaller()



setMethodS3("getAsteriskTags", "PairedPscbsCaller", function(this, collapse=NULL, ...) {
  calls <- this$.calls;
  tags <- c("call", calls);

  # Collapsed or split?
  tags <- Arguments$getTags(tags, collapse=collapse);

  tags;
}, protected=TRUE)


setMethodS3("getRootPath", "PairedPscbsCaller", function(this, ...) {
  "pscbsData";
}, protected=TRUE)


setMethodS3("getPath", "PairedPscbsCaller", function(this, create=TRUE, ...) {
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


setMethodS3("getOptionalArguments", "PairedPscbsCaller", function(this, ...) {
  this$.optionalArgs;
}, protected=TRUE)


setMethodS3("getParameters", "PairedPscbsCaller", function(this, ...) {
  params <- NextMethod("getParameters");
  params$calls <- this$.calls;
  params;
}, protected=TRUE)



setMethodS3("process", "PairedPscbsCaller", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  sms <- getInputDataSet(this);

  verbose && enter(verbose, "Calling LOH and AB");

  if (!force && isDone(this)) {
    verbose && cat(verbose, "Already done. Skipping");
    res <- getOutputDataSet(this);
    verbose && exit(verbose);
    return(res);
  }

  verbose && cat(verbose, "Input data set:");
  verbose && print(verbose, sms);

  pathD <- getPath(this);
  verbose && cat(verbose, "Output path: ", pathD);

  verbose && cat(verbose, "Number of samples: ", length(sms));

  optArgs <- getOptionalArguments(this);
  verbose && cat(verbose, "Optional arguments (may be ignored/may give an error/warning):");
  verbose && str(verbose, optArgs);

  for (ii in seq_along(sms)) {
    smf <- getFile(sms, ii);
    sampleName <- getName(smf);
    verbose && enter(verbose, sprintf("Tumor-normal pair #%d ('%s') of %d", ii, sampleName, length(sms)));

    filename <- getFilename(smf);
    pathname <- file.path(pathD, filename);

    # Sanity check
    stopifnot(getAbsolutePath(pathname) != getAbsolutePath(getFullName(smf)));

    if (!force && isFile(pathname)) {
      verbose && cat(verbose, "Already called. Skipping.");
      verbose && exit(verbose);
      next;
    }

    verbose && enter(verbose, "Loading segmentation data");
    fit <- loadObject(getPathname(smf));

    # Sanity check
    fit <- Arguments$getInstanceOf(fit, "PairedPSCBS");
    verbose && exit(verbose);

    # Arguments to be passed to each caller
    argsT <- append(optArgs, list(verbose=less(verbose, 5)));

    verbose && enter(verbose, "Calling ROH");
    args <- append(list(fit), argsT);
    fit <- do.call(callROH, args);
    verbose && exit(verbose);

    verbose && enter(verbose, "Calling AB");
    args <- append(list(fit), argsT);
    fit <- do.call(callAB, args);
    verbose && exit(verbose);

    verbose && enter(verbose, "Calling LOH");
    args <- append(list(fit), argsT);
    fit <- do.call(callLOH, args);
    verbose && exit(verbose);

    verbose && enter(verbose, "Saving");
    saveObject(fit, file=pathname);
    verbose && exit(verbose);

    verbose && exit(verbose);
  } # for (ii ...)

  res <- getOutputDataSet(this);
  verbose && print(verbose, res);
  verbose && exit(verbose);

  res;
}) # process()


# AD HOC
setMethodS3("getPlatform", "PairedPscbsCaller", function(this, ...) {
  "GenericPlatform";
}, protected=TRUE)



##########################################################################
# HISTORY:
# 2013-01-07
# o BUG FIX: process() for PairedPscbsCaller used the global 'verbose'.
# 2012-09-20
# o Now PairedPscbsCaller() passes '...' to the internal callers.
# 2012-09-19
# o Made as an AromaTransform for now.
# o Created.
##########################################################################
