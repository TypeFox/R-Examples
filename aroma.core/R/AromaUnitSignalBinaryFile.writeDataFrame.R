###########################################################################/**
# @set "class=AromaUnitSignalBinaryFile"
# @RdocMethod writeDataFrame
#
# @title "Writes the data file as a tab-delimited text file"
#
# \description{
#   @get "title" with or without file header comments.
# }
#
# @synopsis
#
# \arguments{
#  \item{filename}{The filename of the generated file.}
#  \item{path}{The path where the generated file should be written.}
#  \item{...}{Not used.}
#  \item{columns}{A @character @vector specifying which column names,
#    including optional annotation data column names, that should be
#    exported.  A string \code{"*"} corresponds to inserting
#    the column names of the source data file.}
#  \item{sep}{A @character string specifying the column separator}.
#  \item{addHeader}{If @TRUE, file header comments will be added that
#    contain additional information about the source file and more.}
#  \item{createdBy}{A @character string specifying the \code{createdBy}
#    file header comment. If @NULL, the package version will be written.}
#  \item{nbrOfDecimals}{An @integer specifying the number of decimals
#    for floating point columns.}
#  \item{sep}{A @character string specifying the column separator}.
#  \item{columnNamesPrefix}{A @character string specifying what type
#    of prefix should be used for column names.}
#  \item{overwrite}{If @TRUE, an existing destination file will be
#    overwritten, otherwise an exception will be thrown.}
#   \item{verbose}{A @logical or a @see "R.utils::Verbose" object.}
# }
#
# \value{
#  Returns the write data file as a @see "TabularTextFile" object.
# }
#
# @author
#
# \seealso{
#  @seeclass
# }
#*/###########################################################################
setMethodS3("writeDataFrame", "AromaUnitSignalBinaryFile", function(this, filename=sprintf("%s.txt", getFilename(this)), path=file.path(getRootName(this, tags="*,txt"), getParentName(this), getChipType(this, fullname=FALSE)), ..., columns=c("unitName", "*"), sep="\t", addHeader=TRUE, createdBy=NULL, nbrOfDecimals=4L, columnNamesPrefix=c("fullname", "name", "none"), overwrite=FALSE, verbose=FALSE) {
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

  fields <- getColumnNames(this);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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

  # Argument 'columnNamesPrefix':
  columnNamesPrefix <- match.arg(columnNamesPrefix);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Writing data frame");

  verbose && cat(verbose, "Data file to be exported:");
  verbose && print(verbose, this);

  verbose && cat(verbose, "Destination pathname: ", pathname);

  # Overwrite?
  if (overwrite && isFile(pathname)) {
    # TODO: Added a backup/restore feature in case new writing fails.
    file.remove(pathname);
    verbose && cat(verbose, "Removed pre-existing file (overwrite=TRUE).");
  }

  # Write to a temporary file
  pathnameT <- pushTemporaryFile(pathname, verbose=verbose);

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


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Generate file header
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfColumns <- nbrOfColumns(this);

  hdr <- NULL;
  if (addHeader) {
    verbose && enter(verbose, "Generating file header");
    hdr <- c(hdr, sprintf("name: %s", getName(this)));
    hdr <- c(hdr, sprintf("tags: %s", getTags(this, collapse=",")));
    hdr <- c(hdr, sprintf("fullName: %s", getFullName(this)));
    hdr <- c(hdr, sprintf("sourceFile: %s", getFilename(this)));
    hdr <- c(hdr, sprintf("sourceChecksum: %s", getChecksum(this)));
    createdOn <- format(Sys.time(), "%Y-%m-%d %H:%M:%S", usetz=TRUE);
    hdr <- c(hdr, sprintf("createdOn: %s", createdOn));
    hdr <- c(hdr, sprintf("createdBy: %s", createdBy));

    hdr <- c(hdr, sprintf("platform: %s", getPlatform(this)));
    hdr <- c(hdr, sprintf("chipType: %s", getChipType(this)));
    hdr <- c(hdr, sprintf("nbrOfUnits: %s", nbrOfUnits(this)));

    # Column types
    columnTypes <- getColumnTypes(this);
    # Sanity check
    columnTypes <- Arguments$getCharacters(columnTypes,
                                      length=rep(nbrOfColumns, times=2L));

    # Annotation data column types
    adColumnTypes <- unname(sapply(adData, FUN=mode));

    # All column types
    allColumnTypes <- c(adColumnTypes, columnTypes);
    allColumnTypesStr <- paste(allColumnTypes, collapse="\t");
    hdr <- c(hdr, sprintf("columnTypes: %s", allColumnTypesStr));


    # Data column names
    columnNames <- getColumnNames(this);
    prefixes <- NULL;
    if (columnNamesPrefix == "name") {
      prefixes <- rep(getName(this), times=nbrOfColumns);
    } else if (columnNamesPrefix == "fullname") {
      prefixes <- rep(getFullName(this), times=nbrOfColumns);
    }
    if (is.null(prefixes)) {
      columnNames <- fields;
    } else {
      if (nbrOfColumns == 1 && fields[1] == "1") {
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
  # Generate data table
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Generating data table");

  # Read data
  verbose && enter(verbose, "Extracting signals");
  data <- this[,,drop=FALSE];
  verbose && exit(verbose);

  # For each data column
  for (cc in seq_len(nbrOfColumns)) {
    values <- data[,cc, drop=TRUE];

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

    data[,cc] <- values;
    # Not needed anymore
    values <- NULL;
  } # for (cc ...)

  verbose && cat(verbose, "Data without annotation data:");
  verbose && str(verbose, data);
  verbose && exit(verbose);

  if (!is.null(adData)) {
    verbose && enter(verbose, "Insert annotation data");
    # Sanity check
    stopifnot(nrow(data) == nrow(adData));
    data <- cbind(adData, data);
    verbose && cat(verbose, "Data with annotation data:");
    verbose && str(verbose, data);
    verbose && exit(verbose);
  }
  # Not needed anymore
  adData <- NULL;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Write to file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Writing file");

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

  write.table(file=pathnameT, data, sep=sep, quote=FALSE, row.names=FALSE,
                                            col.names=FALSE, append=TRUE);
  # Not needed anymore
  data <- NULL;
  verbose && exit(verbose);

  verbose && exit(verbose);


  # Renaming temporary file
  pathname <- popTemporaryFile(pathnameT, verbose=verbose);

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
# 2014-08-27
# o BUG FIX: writeDataFrame(..., columnNamesPrefix="none") for
#   AromaUnitSignalBinaryFile would give an error.
# 2010-05-12
# o Now argument 'path' defaults to <rootPath>,txt/<dataSet>/<chipType>/.
# 2010-04-27
# o Added argument 'columnNamesPrefix'.
# o Removed explicit dependency on aroma.affymetrix. This is now hidden
#   inside getUnitAnnotationDataFile() of AromaPlatformInterface.
# 2010-04-20
# o Removed dependency on getCdf() for affymetrix platforms.
# 2010-04-19
# o Created from aroma.cn::testScripts/system/D9.export/1.export.Rex.
############################################################################
