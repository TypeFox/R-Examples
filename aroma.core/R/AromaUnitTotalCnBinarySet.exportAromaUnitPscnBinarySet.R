###########################################################################/**
# @RdocGeneric exportAromaUnitPscnBinarySet
# @alias exportAromaUnitPscnBinarySet.AromaUnitTotalCnBinarySet
# @alias exportAromaUnitPscnBinarySet.list
#
# @title "Export total and allele B signal data sets as a unified parent-specific copy number signal data set"
#
# \description{
#  @get "title", where each sample is stored in one data file (contrary to the input data sets where each sample is stored in two separated files).
# }
#
# \usage{
#  @usage exportAromaUnitPscnBinarySet,AromaUnitTotalCnBinarySet
#  @usage exportAromaUnitPscnBinarySet,list
# }
#
# \arguments{
#   \item{dsT, dsB}{An @see AromaUnitTotalCnBinarySet and an @see AromaUnitFracBCnBinarySet with coupled sets of samples that match up by name.  If they don't match up, an exception is thrown.  The allele B fraction (BAF) data set \code{dsB} is by default inferred from the total CN data set \code{dsT}.}
#   \item{dataSet, tags}{The name and tags of the output data set.
#     The default is to infer those from the input \code{dsT} data set.}
#   \item{...}{Not used.}
#   \item{rootPath}{The root path of the output data set.}
#   \item{overwrite, skip}{Specifies whether to overwrite and/or skip already exported samples.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns an @see AromaUnitPscnBinarySet object.
# }
#
# \section{Allele-specific CRMAv2 pipeline}{
#   A common use case is to run allele-specific CRMAv2, e.g.
#   \code{dsNList <- doASCRMAv2(csR)}, which outputs a @list \code{dsNList}
#   with elements corresponding to \code{dsT} and \code{dsB}.  This output
#   can be exported to @see AromaUnitPscnBinarySet by this method as
#   \code{dsN <- exportAromaUnitPscnBinarySet(dsNList)}.
# }
#
# @author
#
# @keyword internal
#*/###########################################################################
setMethodS3("exportAromaUnitPscnBinarySet", "AromaUnitTotalCnBinarySet", function(dsT, dsB="*", dataSet="*", tags="*", ..., rootPath="totalAndFracBData/", overwrite=!skip, skip=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dsT':
  dsT <- Arguments$getInstanceOf(dsT, "AromaUnitTotalCnBinarySet");

  # Argument 'dsB':
  if (!identical(dsB, "*")) {
    dsB <- Arguments$getInstanceOf(dsB, "AromaUnitFracBCnBinarySet");
    # Sanity check
    stopifnot(length(dsB) == length(dsT));
    # Reorder
    dsB <- extract(dsB, indexOf(dsB, getNames(dsT)), onDuplicates="error");
    # Sanity check
    stopifnot(getNames(dsB) ==  getNames(dsT));
  }

  # Argument 'tags':
  tags <- Arguments$getTags(tags, collapse=NULL);

  # Argument 'rootPath':
  rootPath <- Arguments$getWritablePath(rootPath);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  clazz <- AromaUnitPscnBinaryFile;

  verbose && enter(verbose, "Exporting to AromaUnitPscnBinarySet");
  if (identical(dsB, "*")) {
    verbose && enter(verbose, "Locating AromaUnitFracBCnBinarySet");
    path <- getPath(dsT);
    verbose && cat(verbose, "Path: ", path);
    pattern <- ",fracB[.]asb$";
    dsB <- AromaUnitFracBCnBinarySet$byPath(path, pattern=pattern);
    verbose && cat(verbose, "Number of files: ", length(dsB));
    dsB <- extract(dsB, indexOf(dsB, getNames(dsT)), onDuplicates="error");
    verbose && cat(verbose, "Number of files kept: ", length(dsB));
    # Sanity check
    stopifnot(length(dsB) == length(dsT));
    verbose && exit(verbose);
  }

  verbose && cat(verbose, "Root path: ", rootPath);
  if (any(dataSet == "*")) {
    default <- getName(dsT, collapse=",");
    dataSet[dataSet == "*"] <- default;
  }
  if (any(tags == "*")) {
    default <- getTags(dsT, collapse=",");
    tags[tags == "*"] <- default;
  }
  dataSetF <- fullname(dataSet, tags=tags);
  verbose && cat(verbose, "Output data set: ", dataSetF);

  chipType <- getChipType(dsT, fullname=FALSE);
  verbose && cat(verbose, "Output chip type: ", chipType);

  path <- filePath(rootPath, dataSetF, chipType);
  path <- Arguments$getWritablePath(path);
  verbose && cat(verbose, "Output path: ", path);

  nbrOfArrays <- length(dsT);
  verbose && cat(verbose, "Number of arrays: ", nbrOfArrays);
  ugp <- getAromaUgpFile(dsT);
  verbose && cat(verbose, "Chip type: ", getChipType(ugp));
  nbrOfUnits <- nbrOfUnits(ugp);

  for (ii in seq_len(nbrOfArrays)) {
    dfT <- dsT[[ii]];
    dfB <- dsB[[ii]];
    name <- getName(dfT);
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d", ii, name, nbrOfArrays));
    # Sanity check
    stopifnot(getName(dfB) == name);
    stopifnot(nbrOfUnits(dfT) == nbrOfUnits);
    stopifnot(nbrOfUnits(dfB) == nbrOfUnits);

    # Convert to original scale
    ratioTag <- NULL;
    logBase <- NULL;
    if (hasTag(dfT, "log2ratio")) {
      logBase <- 2;
      ratioTag <- "ratio";
    } else if (hasTag(dfT, "log10ratio")) {
      logBase <- 10;
      ratioTag <- "ratio";
    } else if (hasTag(dfT, "logRatio")) {
      logBase <- 10;
      ratioTag <- "ratio";
    }

    tagsT <- intersect(getTags(dfT), getTags(dfB));
    tagsT <- c(tagsT, ratioTag);
    fullname <- fullname(name, tags=tagsT);
    filename <- sprintf("%s,pscn.asb", fullname);
    pathname <- Arguments$getWritablePathname(filename, path=path,
                                           mustNotExist=(!overwrite && !skip));
    verbose && cat(verbose, "Pathname: ", pathname);

    if (isFile(pathname)) {
      if (skip) {
        df <- newInstance(clazz, pathname);
        # TODO: We might retrieve an incompatible file.  Validate!
        verbose && cat(verbose, "Already exported. Skipping.");
        verbose && exit(verbose);
        next;
      } else if (!overwrite) {
        throw("Cannot allocate/create file. File already exists: ", pathname);
      }
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Create empty temporary file
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Overwrite?
    if (overwrite && isFile(pathname)) {
      # TODO: Added a backup/restore feature in case new writing fails.
      file.remove(pathname);
      verbose && cat(verbose, "Removed pre-existing file (overwrite=TRUE).");
    }

    pathnameT <- pushTemporaryFile(pathname, verbose=verbose);

    # Read (total,fracB) data
    verbose && enter(verbose, "Reading (total,fracB) data");
    dataT <- extractMatrix(dfT, drop=TRUE);
    dataB <- extractMatrix(dfB, drop=TRUE);
    if (!is.null(logBase)) {
      dataT <- logBase^dataT;
    }
    verbose && exit(verbose);

    # Allocate
    df <- clazz$allocateFromUnitAnnotationDataFile(ugp, filename=pathnameT);

    # Populate
    df[,1L] <- dataT;
    df[,2L] <- dataB;

    # Not needed anymore
    dataT <- dataB <- NULL;

    # Write footer
##    writeFooter(df, footer);

    # Rename temporary file
    pathname <- popTemporaryFile(pathnameT, verbose=verbose);

    # Object to be returned
    df <- newInstance(clazz, pathname);

    verbose && exit(verbose);
  } # for (ii ...)

  verbose && cat(verbose, "Exported data set:");
  pattern <- ",pscn[.]asb$";
  ds <- AromaUnitPscnBinarySet$byPath(path, pattern=pattern);
  ds <- extract(ds, getNames(dsT), onDuplicates="error");
  verbose && print(verbose, ds);

  # Sanity check
  stopifnot(length(ds) == length(dsT));
  stopifnot(all(getNames(ds) == getNames(dsT)));

  verbose && exit(verbose);

  ds;
})


setMethodS3("exportAromaUnitPscnBinarySet", "list", function(dsList, ...) {
  # Argument 'dsList':
  if (length(dsList) != 2L) {
    throw("Argument 'dsList' does not have two elements: ", length(dsList));
  }
  exportAromaUnitPscnBinarySet(dsT=dsList[[1]], dsB=dsList[[2]], ...);
})


############################################################################
# HISTORY:
# 2013-12-11
# o DOCUMENTATION: Added documentation for exportAromaUnitPscnBinarySet().
# 2013-08-12
# o BUG FIX: exportAromaUnitPscnBinarySet() for AromaUnitTotalCnBinarySet
#   would throw an error on "unknown argument 'names' to indexOf()".
# 2012-09-14
# o BUG FIX: exportAromaUnitPscnBinarySet() for AromaUnitTotalCnBinarySet
#   would throw an error if the files of the exported data set was
#   ordered in a non-lexicographic order.
# o Added exportAromaUnitPscnBinarySet() for list.
# 2012-07-21
# o Added exportAromaUnitPscnBinarySet() for AromaUnitTotalCnBinarySet.
# o Created.
############################################################################
