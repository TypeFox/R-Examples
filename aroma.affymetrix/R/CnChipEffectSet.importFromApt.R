###########################################################################/**
# @set "class=CnChipEffectSet"
# @RdocMethod importFromApt
#
# @title "Imports copy-number estimates from an APT summary file"
#
# \description{
#  @get "title".
#  Currently only total copy-number estimates can be imported.
# }
#
# @synopsis
#
# \arguments{
#   \item{filename}{The filename of the APT summary file.}
#   \item{path}{An optional path to the file.}
#   \item{combineAlleles}{If @TRUE, ...}
#   \item{cdf}{An @see "AffymetrixCdfFile" object.}
#   \item{...}{Not used.}
#   \item{skip}{If @TRUE, already imported arrays will be skipped.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns an @see "CnChipEffectSet".
# }
#
# \details{
#   This import method is robust and memory efficient.  One array at the
#   time is imported by first writing to a temporary file which is then
#   renamed to the final name, if import was successful.  (If the import
#   failed, a temporary file will rename that has to be deleted manually).
#
#   Since only one array at the time is imported, the memory overhead
#   will be bounded allowing to import very large tab-delimited data files
#   containing a large number of arrays.  Unfortunately, this approach
#   slows down the reading substantially, because in each import all but
#   one column is parsed but ignored.
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("importFromApt", "CnChipEffectSet", function(static, filename, path=NULL, combineAlleles=TRUE, cdf, ..., skip=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'filename' and 'path':
  pathname <- Arguments$getReadablePathname(filename, path=path, mustExist=TRUE);

  # Argument 'cdf':
  cdf <- Arguments$getInstanceOf(cdf, "AffymetrixCdfFile");

  # Argument 'combineAlleles':
  combineAlleles <- Arguments$getLogical(combineAlleles);
  if (!combineAlleles)
    throw("Currently only 'combineAlleles=TRUE' is supported");

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  ces <- NULL;

  verbose && enter(verbose, "Importing theta estimates from APT summary file");
  verbose && cat(verbose, "Pathname: ", pathname);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Infer data path
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  rootPath <- "plmData";
  chipType <- getChipType(cdf, fullname=FALSE);
  currPath <- dirname(pathname);
  while(TRUE) {
    dataSetName <- basename(currPath);
    if (regexpr(chipType, dataSetName) == -1)
      break;
    currPath <- dirname(currPath);
  }

  outPath <- file.path(rootPath, dataSetName, chipType);
  outPath <- Arguments$getWritablePath(outPath);

  verbose && cat(verbose, "Output path: ", outPath);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get CDF for chip effects
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving CDF for chip effects");
  verbose && printf(verbose, "Chip type: %s,monocell\n", chipType);
   # Get the ChipEffectFile class specific for this set
  clazz <- getChipEffectFileClass(static);
  monocellCdf <- clazz$createParamCdf(cdf);
  # Not needed anymore
  cdf <- NULL;
  verbose && print(verbose, monocellCdf);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Open file & assert file format
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  apt <- AffymetrixAptSummaryFile(pathname);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read header
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  header <- getHeader(apt);
  verbose && cat(verbose, "File header:");
  verbose && str(verbose, header);

  nbrOfArrays <- length(apt);
  arrayNames <- getArrayNames(apt);
  verbose && printf(verbose, "Arrays [%d]: %s\n", nbrOfArrays,
                                        paste(arrayNames, collapse=", "));


  # Infer unit scale
  qScale <- getQuantificationScale(apt);
  if (qScale == "log2") {
  } else {
    # Assume default is log2.
    qScale <- "log2";
  }

  # Garbage collect
  gc <- gc();

  verbose && enter(verbose, "Importing ", nbrOfArrays, " arrays");
  cells <- NULL;
  for (kk in seq_len(nbrOfArrays)) {
    arrayName <- arrayNames[kk];
    verbose && enter(verbose, sprintf("Array #%d (%s) of %d",
                                            kk, arrayName, nbrOfArrays));

    # Create output filename
    filename <- sprintf("%s,chipEffects.CEL", arrayName);
    pathname <- file.path(outPath, filename);

    # Rename lower-case *.cel to *.CEL, if that is the case.  Old versions
    # of the package generated lower-case CEL files. /HB 2007-08-09
    pathname <- AffymetrixFile$renameToUpperCaseExt(pathname);

    verbose && cat(verbose, "Output pathname: ", pathname);

    if (skip && isFile(pathname)) {
      verbose && cat(verbose, "Already imported");
      verbose && exit(verbose);
      next;
    }

    verbose && enter(verbose, "Retrieving chip-effect CEL file");
    verbose && cat(verbose, "Class: ", getName(clazz));

    # Write to a temporary file
    pathnameT <- pushTemporaryFile(pathname, verbose=verbose);

    if (isFile(pathname)) {
      # Test if file can be loaded
      cef <- clazz$fromFile(pathname, verbose=less(verbose));
      # Ok, then rename this file
      res <- file.rename(pathname, pathnameT);
      if (!res) {
        throw("Failed to rename existing file: ", pathname,
                                               " -> ", pathnameT);
      }
      cef <- clazz$fromFile(pathnameT, verbose=less(verbose));
    } else {
      tmpFilename <- basename(pathnameT);
      cef <- clazz$fromDataFile(filename=tmpFilename, path=outPath,
               name=arrayName, cdf=monocellCdf, verbose=less(verbose));
      # Not needed anymore
      tmpFilename <- NULL;
    }
    pathnameT <- getPathname(cef);
    cef$combineAlleles <- combineAlleles;
    cef$mergeStrands <- TRUE;
    verbose && print(verbose, cef);
    verbose && exit(verbose);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);

    # Cell indices in chip-effect CEL files
    if (is.null(cells)) {
      verbose && enter(verbose, "Infering cell indices");

      verbose && enter(verbose, "Reading units");
      dirs <- c("aroma.affymetrix", "APT");
      srcPathname <- getPathname(apt);
      colNames <- getColumnNames(apt);
      key <- list(method="importFromApt", class="CnChipEffectSet", result="importUnits", fileHeader=header, fileSize=file.info(srcPathname)$size, colNames=colNames, chipType=getChipType(monocellCdf));
      verbose && str(verbose, key);

      force <- FALSE;
      res <- NULL;
      if (!force) {
        res <- loadCache(key, dirs=dirs);
        verbose && cat(verbose, "Found cached results");
      }
      if (is.null(res)) {
        # Read unit names
        verbose && enter(verbose, "Retrieving unit names and indices");
        unitNames <- getProbesetIds(apt, verbose=less(verbose, 50));
        # Remove suffix '-A' from SNPs
        isSnp <- grep("SNP_", unitNames);
        unitNames[isSnp] <- gsub("-A$", "", unitNames[isSnp]);
        verbose && str(verbose, unitNames);

        # Match to CDF
        units <- indexOf(cdf, names=unitNames);
        # Not needed anymore
        unitNames <- NULL;
        verbose && str(verbose, units);

        # Any missing?
        keep <- !is.na(units);
        if (any(!keep)) {
          verbose && printf(verbose, "Identified %d unit names not in the CDF.", sum(!keep));
        }
        verbose && exit(verbose);

        res <- list(units=units, keep=keep);
        # Not needed anymore
        units <- keep <- NULL;

        # Store in cache
        saveCache(res, key=key, dirs=dirs);
      }
      units <- res$units;
      keep <- res$keep;
      # Not needed anymore
      res <- NULL;
      # Garbage collect
      gc <- gc();
      verbose && print(verbose, gc);
      verbose && exit(verbose);

      verbose && enter(verbose, "Getting CDF cell indices");
      cells <- lapplyInChunks(units, function(units) {
        cells <- getCellIndices(cef, units=units);
        verbose && cat(verbose, "CDF cell indices:");
        verbose && cat(verbose, "Number of units: ", length(cells));
        cells <- unlist(cells, use.names=FALSE);
        verbose && cat(verbose, "Number of cells: ", length(cells));
        verbose && str(verbose, cells);
        cells <- as.list(cells);
        # Garbage collect
        gc <- gc();
        verbose && print(verbose, gc);
        cells;
      }, chunkSize=100e3, verbose=less(verbose,5));
      cells <- unlist(cells, use.names=FALSE);
      verbose && str(verbose, cells);
      verbose && exit(verbose);

      verbose && exit(verbose);
    } # if (is.null(cells))
    # Not needed anymore
    cef <- NULL;  # Not needed anymore

    verbose && enter(verbose, "Reading data");
    data <- readArrays(apt, arrayName, verbose=less(verbose, 10));
    data <- data[,1];
    if (qScale == "log2") {
      data <- 2^data;
    }
    verbose && str(verbose, data);
    verbose && exit(verbose);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);

    verbose && enter(verbose, "Storing chip effects");
    .updateCel(pathnameT, indices=cells, intensities=data);
    # Not needed anymore
    data <- NULL;
    verbose && exit(verbose);

    # Rename temporary file
    popTemporaryFile(pathnameT, verbose=verbose);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);

    verbose && exit(verbose);
  }
  verbose && exit(verbose);

  # Define chip-effect set
  ces <- byPath(static, path=outPath, cdf=monocellCdf, ...);

  verbose && exit(verbose);

  ces;
}, static=TRUE, private=TRUE)


############################################################################
# HISTORY:
# 2008-01-13
# o Created from importFromDChip().
############################################################################
