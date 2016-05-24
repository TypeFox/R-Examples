###########################################################################/**
# @set "class=AromaUnitSignalBinarySet"
# @RdocMethod writeDataFrame
#
# @title "Writes the data set as a tab-delimited text file"
#
# \description{
#   @get "title" with or without file header comments.
#
#  \emph{We do neither recommend nor encourage the usage of this method;
#   it is only available due to popular demand}.
#   For more details, see below.
# }
#
# @synopsis
#
# \arguments{
#  \item{filename}{The filename of the generated file.}
#  \item{path}{The path where the generated file should be written.}
#  \item{...}{Not used.}
#  \item{units}{The units to be written.  If @NULL, all units are considered.}
#  \item{columns}{A @character @vector specifying which column names,
#    including optional annotation data column names, that should be
#    exported.  A string \code{"*"} corresponds to inserting
#    the column names of the source data files.}
#  \item{columnNamesPrefix}{A @character string specifying what type
#    of prefix should be used for column names.}
#   \item{verbose}{A @logical or a @see "R.utils::Verbose" object.}
# }
#
# \value{
#  Returns the written data set as a @see "TabularTextFileSet" object.
# }
#
# \section{Warning}{
#  There is no limitation in how big the generated file can be.
#  The bigger the data set is, the greater the file size will be.
#  Because of this, \emph{we do neither recommend nor encourage the usage
#  of this method}.  It is available only due to popular demand.
#  Instead we recommend to write each file to a separate file.
#  See \code{writeDataFrame()} of @see AromaUnitSignalBinaryFile for
#  more information.
# }
#
# @author
#
# \seealso{
#  To write the data of one file, see \code{writeDataFrame()} for
#  @see "AromaUnitSignalBinaryFile".
#  @seeclass
# }
#*/###########################################################################
setMethodS3("writeDataFrame", "AromaUnitSignalBinarySet", function(this, filename=sprintf("%s.txt", getFullName(this)), path=file.path(getRootName(this, tags="*,txt"), getFullName(this), getChipType(this, fullname=FALSE)), ..., units=NULL, columns=c("unitName", "*"), sep="\t", addHeader=TRUE, createdBy=NULL, nbrOfDecimals=4L, ram=1, columnNamesPrefix=c("fullnames", "names", "none"), overwrite=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  getAnnotationData <- function(df, columns, ...) {
    res <- list();

    verbose && enter(verbose, "Reading annotation data");

    # Unit genome positions
    verbose && enter(verbose, "Unit genome positions");
    ugp <- getAromaUgpFile(df);
    verbose && cat(verbose, "Unit genome position (UGP) file:");
    verbose && print(verbose, ugp);
    verbose && exit(verbose);
    nbrOfUnits <- nbrOfUnits(ugp);

    # Add unit names?
    if (is.element("unitName", columns)) {
      verbose && enter(verbose, "Unit names");
      unf <- getUnitNamesFile(ugp);
      verbose && cat(verbose, "Filename: ", getFilename(unf));

      unitNames <- getUnitNames(unf);
      verbose && str(verbose, unitNames);
      # Sanity check
      stopifnot(length(unitNames) == nbrOfUnits);

      res$unitName <- unitNames;
      # Not needed anymore
      unitNames <- NULL;

      verbose && exit(verbose);
    }

    if (is.element("chromosome", columns)) {
      res$chromosome <- ugp[,1, drop=TRUE];
    }

    if (is.element("position", columns)) {
      res$position <- ugp[,2, drop=TRUE];
    }

    # Coerce to data frame
    res <- as.data.frame(res, stringsAsFactors=FALSE);

    verbose && cat(verbose, "Annotation data:");
    verbose && str(verbose, res);

    verbose && exit(verbose);

    res;
  } # getAnnotationData()

  knownAnnotationColumnNames <- c("unitName", "chromosome", "position");


  getColumnTypes <- function(this, ...) {
    # Infer column types and column names
    header <- readHeader(this);
    columnTypes <- header$dataHeader$types;
    columnTypes;
  } # getColumnTypes()


  dfFirst <- getOneFile(this, mustExist=TRUE);
  fields <- getColumnNames(dfFirst);
  ugp <- getAromaUgpFile(dfFirst);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'units':
  maxNbrOfUnits <- nbrOfUnits(ugp);
  if (is.null(units)) {
    units <- seq_len(maxNbrOfUnits);
  } else {
    units <- Arguments$getIndices(units, max=maxNbrOfUnits);
  }

  # Argument 'columns':
  columns <- Arguments$getCharacters(columns);

  # Expand asterisk columns
  idxs <- which(columns == "*");
  if (length(idxs) > 0) {
    columns[idxs] <- paste(fields, collapse="\t");
    columns <- strsplit(columns, split="\t", fixed=TRUE);
    columns <- unlist(columns, use.names=FALSE);
  }

  # Check for unknown column names
  knownColumnNames <- c(knownAnnotationColumnNames, fields);
  unknown <- setdiff(columns, knownColumnNames);
  if (length(unknown) > 0) {
    unknown <- paste(unknown, collapse=", ");
    throw("Argument 'columns' contains unknown column names: ", unknown);
  }

  # Argument 'overwrite':
  overwrite <- Arguments$getLogical(overwrite);

  # Argument 'filename' & 'path':
  pathname <- Arguments$getWritablePathname(filename, path=path,
                                                  mustNotExist=!overwrite);

  # Argument 'addHeader':
  addHeader <- Arguments$getLogical(addHeader);

  # Argument 'createdBy':
  if (is.null(createdBy)) {
    pkg <- Package("aroma.core");
    createdBy <- sprintf("%s v%s", getName(pkg), getVersion(pkg));
  } else {
    createdBy <- Arguments$getCharacter(createdBy);
  }

  # Argument 'sep':
  sep <- Arguments$getCharacter(sep);

  # Argument 'nbrOfDecimals':
  nbrOfDecimals <- Arguments$getInteger(nbrOfDecimals, range=c(0,32));

  # Argument 'ram':
  ram <- Arguments$getDouble(ram, range=c(1e-3,999));

  # Argument 'columnNamesPrefix':
  columnNamesPrefix <- match.arg(columnNamesPrefix);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Writing data frame");

  verbose && cat(verbose, "Data set to be exported:");
  verbose && print(verbose, this);

  verbose && cat(verbose, "Destination pathname: ", pathname);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get annotation data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Add annotation data columns?
  adColumnNames <- intersect(knownAnnotationColumnNames, columns);
  if (length(adColumnNames) > 0) {
    adData <- getAnnotationData(this, columns=adColumnNames,
                                                 verbose=less(verbose, 5));
  } else {
    adData <- NULL;
  }
  # Not needed anymore
  adColumnNames <- NULL;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Generate file header
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfFiles <- length(this);
  nbrOfUnits <- length(units);
  nbrOfColumnsPerFile <- nbrOfColumns(dfFirst);
  nbrOfColumns <- nbrOfFiles*nbrOfColumnsPerFile;

  hdr <- NULL;
  if (addHeader) {
    verbose && enter(verbose, "Generating file header");
    hdr <- c(hdr, sprintf("name: %s", getName(this)));
    hdr <- c(hdr, sprintf("tags: %s", getTags(this, collapse=",")));
    hdr <- c(hdr, sprintf("fullName: %s", getFullName(this)));

    hdr <- c(hdr, sprintf("nbrOfFiles: %s", nbrOfFiles));
    hdr <- c(hdr, sprintf("srcFiles: %s", paste(sapply(this, getFilename), collapse="\t")));
    hdr <- c(hdr, sprintf("srcChecksums: %s", paste(sapply(this, getChecksum), collapse="\t")));
    hdr <- c(hdr, sprintf("names: %s", paste(getNames(this), collapse="\t")));
    hdr <- c(hdr, sprintf("tags: %s", paste(sapply(this, getTags, collapse=","), collapse="\t")));
    hdr <- c(hdr, sprintf("fullNames: %s", paste(getFullNames(this), collapse="\t")));

    hdr <- c(hdr, sprintf("platform: %s", getPlatform(this)));
    hdr <- c(hdr, sprintf("chipType: %s", getChipType(this)));
    hdr <- c(hdr, sprintf("nbrOfUnits: %s", nbrOfUnits));

    createdOn <- format(Sys.time(), "%Y-%m-%d %H:%M:%S", usetz=TRUE);
    hdr <- c(hdr, sprintf("createdOn: %s", createdOn));
    hdr <- c(hdr, sprintf("createdBy: %s", createdBy));


    # Data column types
    columnTypes <- getColumnTypes(dfFirst);
    columnTypes <- rep(columnTypes, times=nbrOfFiles);

    # Sanity check
    columnTypes <- Arguments$getCharacters(columnTypes,
                                      length=rep(nbrOfColumns, 2));

    # Annotation data column types
    adColumnTypes <- unname(sapply(adData, FUN=mode));

    # All column types
    allColumnTypes <- c(adColumnTypes, columnTypes);
    allColumnTypesStr <- paste(allColumnTypes, collapse="\t");
    hdr <- c(hdr, sprintf("columnTypes: %s", allColumnTypesStr));


    # Data column names
    prefixes <- NULL;
    if (columnNamesPrefix == "names") {
      prefixes <- rep(getNames(this), each=nbrOfColumnsPerFile);
    } else if (columnNamesPrefix == "fullnames") {
      prefixes <- rep(getFullNames(this), each=nbrOfColumnsPerFile);
    }

    if (is.null(prefixes)) {
      columnNames <- rep(fields, times=nbrOfFiles);
    } else {
      if (nbrOfColumnsPerFile == 1 && fields[1] == "1") {
        columnNames <- prefixes;
      } else {
        columnNames <- paste(prefixes, columnNames, sep=",");
      }
    }

    # Sanity check
    columnNames <- Arguments$getCharacters(columnNames,
                                      length=rep(nbrOfColumns, 2));

    # Annotation data column names
    adColumnNames <- colnames(adData);

    # All column names
    allColumnNames <- c(adColumnNames, columnNames);
    allColumnNamesStr <- paste(allColumnNames, collapse="\t");
    hdr <- c(hdr, sprintf("columnNames: %s", allColumnNamesStr));
    # Not needed anymore
    allColumnNamesStr <- NULL;

    # Turn into file comments
    hdr <- paste("# ", hdr, sep="");
    verbose && cat(verbose, hdr, sep="\n");

    verbose && exit(verbose);
  } # if (addHeader)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Write to file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Writing data set");

  # Overwrite?
  if (overwrite && isFile(pathname)) {
    # TODO: Added a backup/restore feature in case new writing fails.
    file.remove(pathname);
    verbose && cat(verbose, "Removed pre-existing file (overwrite=TRUE).");
  }

  # Write to a temporary file
  pathnameT <- pushTemporaryFile(pathname, verbose=verbose);

  if (!is.null(hdr)) {
    verbose && enter(verbose, "Writing file header");
    cat(file=pathnameT, hdr, sep="\n", append=FALSE);
    verbose && exit(verbose);
  }


  verbose && enter(verbose, "Writing data section");
  quoteColumnNamesIfNeeded <- function(x, quote="\"") {
    idxs <-grep("#", x);
    x[idxs] <- sprintf("%s%s%s", quote, x[idxs], quote);
    x;
  } # quoteColumnNamesIfNeeded()
  allColumnNames <- quoteColumnNamesIfNeeded(allColumnNames);
  allColumnNamesStr <- paste(allColumnNames, collapse=sep);
  cat(file=pathnameT, allColumnNamesStr, sep="\n", append=TRUE);

  # Write in chunks of rows (of constant size)
  chunkSize <- ram * 1e6;
  nbrOfUnitsPerChunk <- chunkSize / nbrOfColumns;
  nbrOfChunks <- ceiling(nbrOfUnits / nbrOfUnitsPerChunk);

  verbose && cat(verbose, "Number of units per chunk: ", nbrOfUnitsPerChunk);
  verbose && cat(verbose, "Number of chunks: ", nbrOfChunks);

  unitsLeft <- units;
  idxsHead <- seq_len(nbrOfUnitsPerChunk);

  verbose && cat(verbose, "Units:");
  verbose && str(verbose, unitsLeft);

  for (kk in seq_len(nbrOfChunks)) {
    verbose && enter(verbose, sprintf("Chunk #%d of %d", kk, nbrOfChunks));

    if (nbrOfUnitsPerChunk > length(unitsLeft)) {
      idxsHead <- seq_len(length(unitsLeft));
    }
    unitsKK <- unitsLeft[idxsHead];

    verbose && cat(verbose, "Units:");
    verbose && str(verbose, unitsKK);

    verbose && enter(verbose, "Extracting signals");
    dataKK <- NULL;
    for (ii in seq_len(nbrOfFiles)) {
      df <- this[[ii]];
      verbose && enter(verbose, sprintf("File #%d ('%s') of %d", ii, getName(df), nbrOfFiles));
      dataII <- df[unitsKK,,drop=FALSE];

      # For each data column
      for (cc in seq_len(nbrOfColumnsPerFile)) {
        values <- dataII[,cc, drop=TRUE];

        # Exporting all missing values as "NA" (not NaN).
        values[is.na(values)] <- NA;
        verbose && cat(verbose, "Signals:");
        verbose && str(verbose, values);

        # Sanity check: Don't allow infinite values
        if (any(is.infinite(values))) {
          throw("Detected infinite values: ", sum(is.infinite(values)));
        }

        # Round off results?
        if (is.double(values) && nbrOfDecimals < Inf) {
          verbose && cat(verbose, "Number of decimals: ", nbrOfDecimals);
          fmtstr <- sprintf("%%.%df", nbrOfDecimals);
          values <- sprintf(fmtstr, values);
        }

        dataII[,cc] <- values;
        # Not needed anymore
        values <- NULL;
      } # for (cc ...)

      if (ii == 1) {
        dataKK <- dataII;
      } else {
        dataKK <- cbind(dataKK, dataII);
      }

      # Not needed anymore
      dataII <- df <- NULL;
      verbose && exit(verbose);
    } # for (ii ...)
    verbose && str(verbose, dataKK);
    verbose && exit(verbose);

    verbose && enter(verbose, "Formatting signals");

    if (!is.null(adData)) {
      verbose && enter(verbose, "Insert annotation data");
      adDataKK <- adData[unitsKK,,drop=FALSE];
      # Sanity check
      stopifnot(nrow(dataKK) == nrow(adDataKK));

      # SUBOPTIMAL: Should preallocate 'dataKK' to save memory
      dataKK <- cbind(adDataKK, dataKK);

      verbose && cat(verbose, "Data with annotation data:");
      verbose && str(verbose, dataKK);
      # Not needed anymore
      adDataKK <- NULL;
      verbose && exit(verbose);
    }

    verbose && str(verbose, dataKK);
    verbose && exit(verbose);

    verbose && enter(verbose, "Writing signals");
    write.table(file=pathnameT, dataKK, sep=sep, quote=FALSE, row.names=FALSE,
                                            col.names=FALSE, append=TRUE);
    # Not needed anymore
    dataKK <- NULL;
    verbose && exit(verbose);

    gc <- gc();
    verbose && print(verbose, gc);

    fileSize <- file.info(pathnameT)$size;
    verbose && printf(verbose, "Current file size: %.1f MB\n", fileSize/1024^2);

    # Next chunk
    unitsLeft <- unitsLeft[-idxsHead];

    verbose && exit(verbose);
  } # for (kk ...)
  # Not needed anymore
  adData <- NULL;

  verbose && exit(verbose);

  # Renaming temporary file
  pathname <- popTemporaryFile(pathnameT, verbose=verbose);

  verbose && exit(verbose);


  verbose && enter(verbose, "Loading output data file");
  res <- TabularTextFile(pathname, sep=sep);
  verbose && print(verbose, res);
  verbose && exit(verbose);

  verbose && exit(verbose);

  # Return results
  res;
}) # writeDataFrame()


############################################################################
# HISTORY:
# 2011-12-15 [HB]
# o Added argument 'units' to writeDateFrame() for
#   AromaUnitSignalBinarySet to make it possible to write any subset
#   of units and in any order (e.g. genome order).
# 2010-07-07 [PN]
# o BUG FIX: writeDateFrame() for AromaUnitSignalBinarySet would write the
#   same data chunk over and over.
# 2010-05-12
# o Now argument 'path' defaults to <rootPath>,txt/<dataSet>/<chipType>/.
# o BUG FIX: writeDateFrame() for AromaUnitSignalBinarySet was not working.
# 2010-04-27
# o Added argument 'columnNamesPrefix'.
# o Removed explicit dependency on aroma.affymetrix. This is now hidden
#   inside getUnitAnnotationDataFile() of AromaPlatformInterface.
# 2010-04-22
# o Created AromaUnitSignalBinaryFile.writeDataFrame.R.
############################################################################
