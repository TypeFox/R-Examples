###########################################################################/**
# @RdocClass AromaTabularBinaryFile
#
# @title "The AromaTabularBinaryFile class"
#
# \description{
#  @classhierarchy
#
#  A AromaTabularBinaryFile represents a file with a binary format.
#  It has a well defined header, a data section, and a footer.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "R.filesets::GenericTabularFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#
# \seealso{
#   @see "R.filesets::GenericDataFile".
# }
#*/###########################################################################
setConstructorS3("AromaTabularBinaryFile", function(...) {
  this <- extend(GenericTabularFile(..., .verify=FALSE),
                   c("AromaTabularBinaryFile", uses("FileCacheKeyInterface")),
    "cached:.hdr"=NULL,
    "cached:.ftr"=NULL
  );

  # Parse attributes (all subclasses must call this in the constructor).
  pathname <- getPathname(this)
  setAttributesByTags(this)

  this;
})


setMethodS3("as.character", "AromaTabularBinaryFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character");
  s <- c(s, sprintf("File format: v%d", readHeader(this)$fileVersion));
  s <- c(s, sprintf("Dimensions: %dx%d", nbrOfRows(this), nbrOfColumns(this)));
  s <- c(s, sprintf("Column classes: %s",
                         paste(getColClasses(this), collapse=", ")));
  s <- c(s, sprintf("Number of bytes per column: %s",
                         paste(getBytesPerColumn(this), collapse=", ")));

  footer <- readFooter(this, asXmlString=TRUE);
  footer <- gsub(">[\n\r ]*", ">", footer);
  footer <- gsub("^[ ]*", "", footer);
  s <- c(s, sprintf("Footer: %s", footer));
  s;
}, protected=TRUE)


setMethodS3("setAttributesByTags", "AromaTabularBinaryFile", function(this, ...) {
  # Does nothing.
}, protected=TRUE)


setMethodS3("getDefaultColumnNames", "AromaTabularBinaryFile", function(this, ...) {
  as.character(seq_len(nbrOfColumns(this)));
}, protected=TRUE)


setMethodS3("dimnames<-", "AromaTabularBinaryFile", function(x, value) {
  # To please R CMD check
  this <- x;

  throw("Dimension names of an ", class(this)[1], " are read only.");
}, appendVarArgs=FALSE, protected=TRUE)


setMethodS3("readHeader", "AromaTabularBinaryFile", function(this, con=NULL, ..., force=FALSE) {
  if (is.null(con)) {
    # Look for cached results
    hdr <- this$.hdr;
    if (!force && !is.null(hdr))
      return(hdr);
  }

  knownDataTypes <- c("integer"=1, "double"=2, "raw"=3);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  readBytes <- function(con, n=1, ...) {
    readBin(con=con, what=integer(), size=1, n=n, signed=FALSE, endian="little");
  }

  readShorts <- function(con, n=1, ...) {
    readBin(con=con, what=integer(), size=2, n=n, signed=FALSE, endian="little");
  }

  # Non-signed integers of length 4 bytes are not supported, cf. help(readBin).
  readInts <- function(con, n=1, ...) {
    readBin(con=con, what=integer(), size=4, n=n, signed=TRUE, endian="little");
  }

  readString <- function(con, ...) {
    nbrOfBytes <- readInts(con);
    nbrOfBytes <- Arguments$getInteger(nbrOfBytes, range=c(0,2^20));
    readChar(con=con, nchars=nbrOfBytes);
  }

  readDataHeader <- function(con, ...) {
    # Number of elements (rows)
    nbrOfRows <- readInts(con);
    nbrOfRows <- Arguments$getInteger(nbrOfRows, range=c(0,200e6));

    # Number of fields (columns)
    nbrOfColumns <- readInts(con);
    nbrOfColumns <- Arguments$getInteger(nbrOfColumns, range=c(0,1000));

    # Types of columns
    types <- readBytes(con, n=nbrOfColumns);
    types <- Arguments$getIntegers(types, range=range(knownDataTypes));
    types <- names(knownDataTypes)[types];

    # Number of bytes per column
    sizes <- readBytes(con, n=nbrOfColumns);
    sizes <- Arguments$getIntegers(sizes, range=c(1,8));
    ok <- (sizes %in% c(1,2,4,8));
    if (any(!ok)) {
      cc <- which(!ok);
      throw("File format error. Detect one or more columns with invalid byte sizes, i.e. not in {1,2,4,8}: ", paste(paste(cc, sizes[cc], sep=":"), collapse=", "));
    }

    # Assert that 'raw' columns are only of size one
    nok <- (sizes[types == "raw"] != 1);
    if (any(nok)) {
      cc <- which(nok);
      throw("File format error. Detect one or more columns of data type 'raw' but of size different from one: ", paste(paste(cc, sizes[cc], sep=":"), collapse=", "));
    }

    # Are the columns signed or not?
    signeds <- readBytes(con, n=nbrOfColumns);
    signeds <- Arguments$getIntegers(signeds, range=c(0,1));
    signeds <- as.logical(signeds);

    nbrOfBytes <- nbrOfRows*sizes;
    dataOffsets <- c(0, cumsum(nbrOfBytes[-length(nbrOfBytes)]));

    dataOffset <- seek(con=con, rw="r");

    # Offset to the footer, which follows immediately after the data
    # section.
    footerOffset <- dataOffset + dataOffsets[nbrOfColumns] +
                                                 nbrOfBytes[nbrOfColumns];
    list(
      nbrOfRows=nbrOfRows,
      nbrOfColumns=nbrOfColumns,
      types=types,
      sizes=sizes,
      signeds=signeds,
      nbrOfBytes=nbrOfBytes,
      dataOffsets=dataOffsets,
      dataOffset=dataOffset,
      footerOffset=footerOffset
    );
  } # readDataHeader()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Main
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Open file?
  if (is.null(con)) {
    pathname <- getPathname(this);
    con <- file(pathname, open="rb");
    on.exit(close(con));
  }

  # Read magic
  trueMagic <- charToRaw("aroma");
  magic <- readBin(con=con, what=raw(), n=length(trueMagic));
  if (!identical(magic, trueMagic)) {
    asStr <- function(raw) {
      paste("[", paste(sprintf("%#0x", as.integer(raw)), collapse=","),
                                                             "]", sep="");
    }
    throw("File format error. The read \"magic\" does not match the existing one: ", asStr(magic), " != ", asStr(trueMagic));
  }

  version <- readInts(con, n=1);
  if (version < 0) {
    throw("File format error. Negative file version: ", version);
  }
  if (version > 10e3) {
    throw("File format error. Ridiculously large file version (>10e3): ", version);
  }

  if (version >= 1 && version <= 1) {
    comment <- readString(con=con);
    dataHeader <- readDataHeader(con=con);
  } else {
    throw("Unknown file format version: ", version);
  }

  hdr <- list(fileVersion=version, comment=comment, dataHeader=dataHeader);

  # Cache result
  this$.hdr <- hdr;

  hdr;
}, protected=TRUE)


setMethodS3("readRawFooter", "AromaTabularBinaryFile", function(this, con=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Non-signed integers of length 4 bytes are not supported, cf. help(readBin)
  readInts <- function(con, n=1, ...) {
    readBin(con=con, what=integer(), size=4, n=n, signed=TRUE, endian="little");
  }

  if (is.null(con)) {
    # Look for cached results
    ftr <- this$.ftr;
    if (!is.null(ftr))
      return(ftr);
  }

  # Open file?
  if (is.null(con)) {
    pathname <- getPathname(this);
    con <- file(pathname, open="rb");
    on.exit(close(con));
  }

  hdr <- readHeader(this, con=con, ...);
  footerOffset <- hdr$dataHeader$footerOffset;

  # Move to the footer
  seek(con=con, where=footerOffset, origin="start", rw="r");
  nbrOfBytes <- readInts(con=con, size=4);

  raw <- readBin(con=con, what=raw(), n=nbrOfBytes);

  res <- list(
    nbrOfBytes=nbrOfBytes,
    raw=raw
  );

  res;
}, protected=TRUE)



###########################################################################/**
# @RdocMethod readFooter
#
# @title "Reads the file footer in XML format into a named nested list"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{asXmlString}{If @TRUE, the file footer is returned as
#      a @character string.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a named @list structure (or a @character string).
# }
#
# @author
#
# \seealso{
#   @seemethod "writeFooter".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("readFooter", "AromaTabularBinaryFile", function(this, asXmlString=FALSE, ...) {
  raw <- readRawFooter(this)$raw;
  if (length(raw) == 0) {
    if (asXmlString)
      return("");
    return(NULL);
  }

  xml <- rawToChar(raw);
  if (asXmlString) {
    xml <- trim(xml);
    xml <- gsub("^<footer>", "", xml);
    xml <- trim(xml);
    xml <- gsub("</footer>$", "", xml);
    xml <- trim(xml);
    res <- xml;
  } else {
    res <- xmlToList(xml);
    if (identical(names(res), "footer"))
      res <- res[["footer"]];
  }
  res;
})





###########################################################################/**
# @RdocMethod writeFooter
#
# @title "Writes a named nested list to the file footer in XML format"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{footer}{A named @list structure.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns nothing.
# }
#
# @author
#
# \seealso{
#   @seemethod "readFooter".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("writeFooter", "AromaTabularBinaryFile", function(this, footer, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'footer':
  if (!is.list(footer)) {
    throw("Argument 'footer' is not a list: ", mode(footer));
  }
  if (identical(names(footer), "footer")) {
    footer <- list(footer=footer);
  }

  # Generate XML version of attributes
  xml <- listToXml(footer, indentStep="");
  xml <- trim(xml);

  # Generate raw byte stream of attributes
  raw <- charToRaw(xml);

  # Write to file
  writeRawFooter(this, raw);
})


setMethodS3("writeRawFooter", "AromaTabularBinaryFile", function(this, raw, con=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Locale functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  writeInts <- function(con, values, ...) {
    values <- as.integer(values);
    writeBin(con=con, values, size=4, endian="little");
  }

  # Non-signed integers of length 4 bytes are not supported, cf. help(readBin)
  readInts <- function(con, n=1, ...) {
    readBin(con=con, what=integer(), size=4, n=n, signed=TRUE, endian="little");
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Writing footer");

  # Need to open a file?
  isFile <- (is.null(con));
  if (isFile) {
    pathname <- getPathname(this);
    # Sanity check
    stopifnot(isFile(pathname));
    pathname <- Arguments$getWritablePathname(pathname);
    con <- file(pathname, open="r+b");
    verbose && cat(verbose, "Opened file ('r+b') to be close automatically");
    verbose && cat(verbose, "Pathname: ", pathname);
    on.exit(close(con), add=TRUE);
  }


  verbose && enter(verbose, "Locating footer");
  hdr <- readHeader(this, con=con, ...);
  footerOffset <- hdr$dataHeader$footerOffset;
  verbose && cat(verbose, "File position: ", footerOffset);

  # Read current footer
  seek(con=con, where=footerOffset, origin="start", rw="r");
  nbrOfBytes <- readInts(con=con, size=4);
  verbose && cat(verbose, "Current length of footer: ", nbrOfBytes);
  verbose && exit(verbose);


  verbose && enter(verbose, "Modifying footer");
  # Move to the footer
  seek(con=con, where=footerOffset, origin="start", rw="w");

  # Write length of footer
  size <- length(raw);
  writeInts(con=con, size);
  writeBin(con=con, raw);
  verbose && enter(verbose, "New length: ", size);

  # Erase the rest of the footer
  rest <- nbrOfBytes - size;
  if (rest > 0) {
    verbose && enter(verbose, "Clearing reminder of footer");
    verbose && cat(verbose, "Number of bytes: ", rest);
    writeBin(con=con, raw(rest));
    verbose && exit(verbose);
  }
  verbose && exit(verbose);

  verbose && exit(verbose);
}, protected=TRUE)


setMethodS3("readDataFrame", "AromaTabularBinaryFile", function(this, rows=NULL, columns=NULL, ..., retRowNames=FALSE, drop=FALSE, verbose=FALSE) {
  # Open file
  pathname <- getPathname(this);
  con <- file(pathname, open="rb");
  on.exit(close(con));

  # Data header
  hdr <- readHeader(this, con=con)$dataHeader;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'rows':
  rownames <- NULL;
  if (is.null(rows)) {
    rows <- seq_len(hdr$nbrOfRows);
  } else if (is.logical(rows)) {
    rows <- which(rows);
    rows <- Arguments$getIndices(rows, max=hdr$nbrOfRows);
    if (retRowNames) {
      rownames <- as.character(rows);
      rownames <- make.unique(rownames);
    }
  } else {
    rows <- Arguments$getIndices(rows, max=hdr$nbrOfRows);
    if (retRowNames) {
      rownames <- as.character(rows);
      rownames <- make.unique(rownames);
    }
  }

  # Argument 'columns':
  if (is.null(columns)) {
    columns <- seq_len(hdr$nbrOfColumns);
  } else if (is.logical(columns)) {
    columns <- which(columns);
    columns <- Arguments$getIndices(columns, max=hdr$nbrOfColumns);
  } else {
    columns <- Arguments$getIndices(columns, max=hdr$nbrOfColumns);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose), add=TRUE);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Allocate return object
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Allocating data object");
  colClasses <- hdr$types[columns];
  verbose && cat(verbose, "Number of rows: ", length(rows));
  verbose && cat(verbose, "Column classes: ", paste(colClasses, collapse=", "));
  data <- dataFrame(colClasses=colClasses, nrow=length(rows));
  if (!is.null(rownames))
    rownames(data) <- rownames;
  verbose && str(verbose, data, level=-30);
  verbose && exit(verbose);

  # Nothing more todo?
  if (length(rows) == 0) {
    colnames(data) <- getColumnNames(this)[columns];

    # Drop singleton dimensions?
    if (drop) {
      if (ncol(data) == 1) {
        data <- data[,1,drop=TRUE];
      }
    }

    return(data);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Reading data");
  # First and last row to read in each column
  rrow <- range(rows);
  nbrOfRows <- as.integer(diff(rrow)+1);
  # Shift rows such that min(rows) == 1.
  rows <- rows - as.integer(rrow[1] - 1);

  # Record the current file offset
  dataOffsets <- hdr$dataOffsets[columns];

  # Read data in the order it appears on file
  o <- order(dataOffsets);

  count <- 0;
  for (kk in o) {
    count <- count + 1;
    verbose && enter(verbose, "Reading column #", count, " of ", length(o), level=-20);
    cc <- columns[kk];
    type <- hdr$types[cc];
    size <- hdr$sizes[cc];
    signed <- hdr$signeds[cc];

    verbose && printf(verbose, "Column %d: %s, %d bytes, signed=%s\n", cc, type, size, signed, level=-50);

    # Jump to the start of the data block
    dataOffset <- hdr$dataOffset + dataOffsets[kk] + (rrow[1]-1)*size;
    verbose && printf(verbose, "Data offset: %d\n", dataOffset, level=-50);
    seek(con=con, where=dataOffset, origin="start", rw="r");

    # Read from first to last row to be read, the discard unwanted.
    # TO DO: Optimize this.
    verbose && enter(verbose, "Reading binary data", level=-20);
    values <- readBin(con=con, n=nbrOfRows, what=type, size=size,
                                         signed=signed, endian="little");
    verbose && exit(verbose);

    values <- values[rows];

    # Store data
    data[[o[kk]]] <- values;

    # Not needed anymore
    values <- NULL;
    verbose && exit(verbose);
  }
  verbose && exit(verbose);

  # Add column names
  colnames(data) <- getColumnNames(this)[columns];

  # Drop singleton dimensions?
  if (drop) {
    if (ncol(data) == 1) {
      data <- data[,1,drop=TRUE];
    } else if (nrow(data) == 1) {
      data <- data[1,,drop=TRUE];
    }
  }

  data;
}, protected=TRUE)



setMethodS3("readColumns", "AromaTabularBinaryFile", function(this, ...) {
  readDataFrame(this, ...);
})


setMethodS3("updateDataColumn", "AromaTabularBinaryFile", function(this, rows=NULL, column, values, .con=NULL, .hdr=NULL, .validateArgs=TRUE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  con <- .con;
  if (!is.null(con))
    seek(con, where=0, offset="start", rw="r");
  hdr <- .hdr;

  if (.validateArgs) {
    if (is.null(con)) {
      # Open file
      pathname <- getPathname(this);
      # Sanity check
      stopifnot(isFile(pathname));
      pathname <- Arguments$getWritablePathname(pathname);
      con <- file(pathname, open="r+b");
      on.exit(close(con));
    }

    # Data header
    if (is.null(hdr)) {
      hdr <- readHeader(this, con=con)$dataHeader;
    }

    # Argument 'rows':
    if (is.null(rows)) {
      rows <- seq_len(hdr$nbrOfRows);
    } else if (is.logical(rows)) {
      rows <- which(rows);
      rows <- Arguments$getIndices(rows, max=hdr$nbrOfRows);
    } else {
      rows <- Arguments$getIndices(rows, max=hdr$nbrOfRows);
    }

    # Argument 'column':
    column <- Arguments$getIndex(column, max=hdr$nbrOfColumns);

    # Argument 'verbose':
    verbose <- Arguments$getVerbose(verbose);
    if (verbose) {
      pushState(verbose);
      on.exit(popState(verbose), add=TRUE);
    }
  } # if (.validateArgs)

  verbose && enter(verbose, "Updating data column by writing to file");

  verbose && cat(verbose, "Number of rows: ", length(rows));
  verbose && cat(verbose, "Column: ", column);
  verbose && printf(verbose, "Values: %d %s(s)\n", length(values), mode(values));

  values <- rep(values, length.out=length(rows));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Prepare data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Optimizing data to be writing");
  verbose && cat(verbose, "Rows and values:");
  verbose && str(verbose, rows);
  verbose && str(verbose, values);

  # Remove duplicated rows
  rows <- rev(rows);
  values <- rev(values);
  dups <- duplicated(rows);
  rows <- rows[!dups];
  values <- values[!dups];
  # Not needed anymore
  dups <- NULL;

  # Reorder rows
  o <- order(rows);
  rows <- rows[o];
  values <- values[o];
  # Not needed anymore
  o <- NULL;

  type <- hdr$types[column];
  size <- hdr$sizes[column];
  signed <- hdr$signeds[column];

  # Censor raw and integer data
  if (type %in% c("raw", "integer")) {
    if (type == "raw") {
      range <- c(0, 255);
    } else if (type == "integer") {
      # FYI: intNA <- as.integer(2^31);
      if (signed) {
        range <- c(-2^(8*size-1), 2^(8*size-1)-1);
      } else {
        range <- c(0, 2^(8*size)-1);
      }
    }

    # Coerce values?
    if (!mode(values) %in% c("raw", "numeric")) {
      values <- as.double(values);
    }

    msgL <- msgH <- NULL;

    idxs <- which(values < range[1]);
    nL <- length(idxs);
    if (nL > 0) {
      rangeL <- range(values[idxs], na.rm=TRUE);
      msgL <- sprintf("%d values in [%.0f,%.0f] were too small",
                                       nL, rangeL[1], rangeL[2]);
      values[idxs] <- range[1];
    }
    idxs <- which(values > range[2]);
    nH <- length(idxs);
    if (nH > 0) {
      rangeH <- range(values[idxs], na.rm=TRUE);
      msgH <- sprintf("%d values in [%.0f,%.0f] were too large",
                                       nH, rangeH[1], rangeH[2]);
      values[idxs] <- range[2];
    }

    if (nL+nH > 0) {
      msg <- sprintf("%d values to be assigned were out of range [%.0f,%.0f] and therefore censored to fit the range. Of these, %s.", (nL+nH), range[1], range[2], paste(c(msgL, msgH), collapse=" and "));
      verbose && cat(verbose, msg);
      warning(msg);
    }
  }


  # Coerce data
  # Data type information
  storage.mode(values) <- type;

  verbose && cat(verbose, "Rows and values:");
  verbose && str(verbose, rows);
  verbose && str(verbose, values);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Write data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Shift rows such that min(rows) == 1.
  firstRow <- rows[1];
  rows <- rows - firstRow + 1;
  nbrOfRows <- rows[length(rows)];

  # Calculate the offset of the first element to read/write
  dataOffset <- hdr$dataOffset + hdr$dataOffsets[column] + (firstRow-1)*size;

  # 1) Read existing data
  verbose && enter(verbose, "Reading existing data");
  seek(con=con, where=dataOffset, origin="start", rw="r");
  signed <- hdr$signeds[column];
  oldValues <- readBin(con=con, n=nbrOfRows, what=type, size=size, signed=signed, endian="little");
  verbose && str(verbose, oldValues);
  verbose && exit(verbose);

  # 2) Coerce and update the values
  storage.mode(oldValues) <- type;
  oldValues[rows] <- values;
  verbose && str(verbose, oldValues);
  # Not needed anymore
  values <- rows <- NULL;

  # 3) Write back
  verbose && enter(verbose, "Writing updated data");
  seek(con=con, where=dataOffset, origin="start", rw="w");
  writeBin(con=con, object=oldValues, size=size, endian="little");
  flush(con);
  verbose && exit(verbose);

  verbose && exit(verbose);

  invisible(this);
}, protected=TRUE)



setMethodS3("updateData", "AromaTabularBinaryFile", function(this, rows=NULL, columns=NULL, values, ..., verbose=FALSE) {
  # Open file
  pathname <- getPathname(this);
  # Sanity check
  stopifnot(isFile(pathname));
  pathname <- Arguments$getWritablePathname(pathname);
  con <- file(pathname, open="r+b");
  on.exit(close(con));

  # Data header
  hdr <- readHeader(this, con=con)$dataHeader;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'rows':
  if (is.null(rows)) {
    rows <- seq_len(hdr$nbrOfRows);
  } else if (is.logical(rows)) {
    rows <- which(rows);
    rows <- Arguments$getIndices(rows, max=hdr$nbrOfRows);
  } else {
    rows <- Arguments$getIndices(rows, max=hdr$nbrOfRows);
  }
  nbrOfRows <- length(rows);

  # Argument 'columns':
  if (is.null(columns)) {
    columns <- seq_len(hdr$nbrOfColumns);
  } else if (is.logical(columns)) {
    columns <- which(columns);
    columns <- Arguments$getIndices(columns, max=hdr$nbrOfColumns);
  } else {
    columns <- Arguments$getIndices(columns, max=hdr$nbrOfColumns);
  }
  nbrOfColumns <- length(columns);

  # Argument 'values':
  if (is.data.frame(values) || is.matrix(values)) {
    if (ncol(values) != nbrOfColumns) {
      throw("Number of columns in ", class(values), " 'values' does not match the number of specified columns: ", ncol(values), " != ", nbrOfColumns);
    }
  } else if (is.list(values)) {
    if (length(values) != nbrOfColumns) {
      throw("Number of elements in list 'values' does not match the number of specified columns: ", length(values), " != ", nbrOfColumns);
    }
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose), add=TRUE);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Update each column
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Update the column in order, because that is faster
  o <- order(columns);
  count <- 0;
  for (kk in o) {
    count <- count + 1;
    verbose && enter(verbose, "Updating column #", count, " of ", length(o));
    cc <- o[kk];
    column <- columns[cc];

    # Extract the values
    if (is.data.frame(values) || is.matrix(values)) {
      theValues <- values[,cc];
    } else if (is.list(values)) {
      theValues <- values[[cc]];
    } else {
      # Is this strange?
      theValues <- values;
    }

    updateDataColumn(this, .con=con, .hdr=hdr, rows=rows, column=column, values=theValues, verbose=less(verbose));

    # Not needed anymore
    theValues <- NULL;

    verbose && exit(verbose);
  } # for (kk ...)

  invisible(this);
}, protected=TRUE)



###########################################################################/**
# @RdocMethod allocate
#
# @title "Creates an AromaTabularBinaryFile"
#
# \description{
#  @get "title" of a certain dimension and data column types.
# }
#
# @synopsis
#
# \arguments{
#   \item{filename}{The filename of the new file.}
#   \item{path}{The path where to store the new file.}
#   \item{nbrOfRows}{An @integer specifying the number of rows to allocate.}
#   \item{types}{A @character @vector specifying the data type of each
#      column.  The length specifies the number of columns to allocate.}
#   \item{sizes}{An @integer @vector of values in \{1,2,4,8\} specifying
#      the size of each column (data type).}
#   \item{signeds}{An @logical @vector specifying if the data types in each
#      column is signed or not.}
#   \item{defaults}{An optional @list (or @vector) containing default
#      values for each column.}
#   \item{comment}{An optional @character string written to the file header.}
#   \item{overwrite}{If @TRUE, an existing file is overwritten, otherwise not.}
#   \item{skip}{If @TRUE and \code{overwrite=TRUE}, any existing file is
#      returned as is.}
#   \item{footer}{An optional @list of attributes written (as character
#      strings) to the file footer.}
#   \item{...}{Not used.}
#   \item{verbose}{@see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @see "AromaTabularBinaryFile" object.
# }
#
# \section{Data types}{
#   Valid data types are currently "@integer" and "@double".
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("allocate", "AromaTabularBinaryFile", function(static, filename, path=NULL, nbrOfRows, types, sizes, signeds=TRUE, defaults=NA, comment=NULL, overwrite=FALSE, skip=FALSE, footer=list(), ..., verbose=FALSE) {
  knownDataTypes <- c("integer"=1, "double"=2, "raw"=3);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  writeBytes <- function(con, values, ...) {
    values <- as.integer(values);
    writeBin(con=con, values, size=1, endian="little");
  }

  writeShorts <- function(con, values, ...) {
    values <- as.integer(values);
    writeBin(con=con, values, size=2, endian="little");
  }

  writeInts <- function(con, values, ...) {
    values <- as.integer(values);
    writeBin(con=con, values, size=4, endian="little");
  }

  writeFloats <- function(con, values, ...) {
    values <- as.double(values);
    writeBin(con=con, values, size=4, endian="little");
  }

  writeDoubles <- function(con, values, ...) {
    values <- as.double(values);
    writeBin(con=con, values, size=8, endian="little");
  }

  writeString <- function(con, value, ...) {
    writeInts(con, nchar(value));  # Note, it is NOT an zero-terminated string
    writeChar(con=con, value, nchars=nchar(value), eos=NULL);
  }


  writeDataHeader <- function(con, nbrOfRows, types, sizes, signeds, ...) {
    # Number of elements (rows)
    writeInts(con=con, nbrOfRows);

    # Number of fields (columns)
    nbrOfColumns <- length(types);
    writeInts(con=con, nbrOfColumns);

    # Types of columns
    types <- knownDataTypes[types];
    writeBytes(con=con, types);

    # Number of bytes per column
    writeBytes(con=con, sizes);

    # Are the columns signed or not?
    writeBytes(con=con, signeds);
  } # writeDataHeader()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'nbrOfRows':
  nbrOfRows <- Arguments$getInteger(nbrOfRows, range=c(0,200e6));

  # Argument 'types':
  if (is.character(types))
    types <- knownDataTypes[types];
  types <- Arguments$getIntegers(types, range=range(knownDataTypes));
  types <- names(knownDataTypes)[types];

  nbrOfColumns <- length(types);
  nbrOfColumns <- Arguments$getInteger(nbrOfColumns, range=c(0,1000));

  # Argument 'sizes':
  sizes <- Arguments$getIntegers(sizes, range=c(1,8));
  ok <- (sizes %in% c(1,2,4,8));
  if (any(!ok)) {
    cc <- which(!ok);
    throw("Cannot allocate/create file. Detect one or more columns with invalid byte sizes, i.e. not in {1,2,4,8}: ", paste(paste(cc, sizes[cc], sep=":"), collapse=", "));
  }
  sizes <- rep(sizes, length.out=nbrOfColumns);

  # Check (types, sizes)
  if (any(types == "raw" & sizes > 1)) {
    throw("Raws can only be stored as single bytes.");
  }
  if (any(types == "integer" & sizes > 4)) {
    throw("Integers can only be stored as 1, 2 or 4 bytes, not 8.");
  }
  if (any(types == "integer" & sizes == 4 & !signeds)) {
    throw("Integers stored in 4 bytes must be signed.");
  }

  # Argument 'signeds':
  signeds <- Arguments$getLogicals(signeds);
  signeds <- rep(signeds, length.out=nbrOfColumns);

  # Argument 'defaults':
  defaults <- rep(defaults, length.out=nbrOfColumns);
  defaults <- as.list(defaults);
  for (kk in seq_along(defaults)) {
    storage.mode(defaults[[kk]]) <- types[kk];
  }

  # Argument 'comment':
  if (is.null(comment)) {
    pkg <- "aroma.core";
    ver <- packageDescription(pkg)$Version;
    comment <- sprintf("Created by the %s (v%s) package.", pkg, ver);
  }

  # Argument 'path':
  path <- Arguments$getWritablePath(path);


  # Argument 'footer':
  if (is.null(footer)) {
  } else if (!is.list(footer)) {
    throw("Argument 'footer' must be NULL or a list: ", class(footer)[1]);
  }

  # Argument 'overwrite':
  overwrite <- Arguments$getLogical(overwrite);

  # Argument 'skip':
  skip <- Arguments$getLogical(skip);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Main
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pathname <- Arguments$getWritablePathname(filename, path=path,
                                         mustNotExist=(!overwrite && !skip));
  verbose && cat(verbose, "Pathname: ", pathname);

  if (isFile(pathname)) {
    if (skip) {
      res <- newInstance(static, pathname);
      # TODO: We might retrieve an incompatible file.  Validate!
      return(res);
    } else if (!overwrite) {
      throw("Cannot allocate/create file.  File already exists: ", pathname);
    }
  }

  verbose && cat(verbose, "nbrOfRows: ", nbrOfRows);
  verbose && cat(verbose, "nbrOfColumns: ", nbrOfColumns);
  verbose && cat(verbose, "types: ", paste(types, collapse=", "));
  verbose && cat(verbose, "sizes: ", paste(sizes, collapse=", "));
  verbose && cat(verbose, "signed: ", paste(signeds, collapse=", "));
  verbose && cat(verbose, "defaults:");
  verbose && str(verbose, defaults);

  verbose && cat(verbose, "Attributes to be written to file footer:");
  verbose && str(verbose, footer);

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

  con <- file(pathnameT, open="wb");
  on.exit({
    if (!is.null(con))
      close(con);
    con <- NULL;
  }, add=TRUE);

  # Write magic
  magic <- charToRaw("aroma");
  writeBin(con=con, magic);

  # Write file version
  version <- 1;
  writeInts(con=con, version);

  # Write comment
  writeString(con=con, comment);

  # Write data header
  writeDataHeader(con=con, nbrOfRows=nbrOfRows, types=types, sizes=sizes, signeds=signeds);

  # Write empty data, column by column
  for (cc in seq_len(nbrOfColumns)) {
    size <- sizes[cc];
    type <- types[cc];
    signed <- signeds[cc];
    default <- defaults[[cc]];
    values <- rep(default, times=nbrOfRows);
    writeBin(con=con, values, size=size, endian="little");
    # Not needed anymore
    values <- NULL;
  }

  # Write empty footer (this may be used to add extra meta data)
  # Write size of footer
  size <- 0;
  writeInts(con=con, size);

  # Close connection (otherwise writeFooter() will fail below)
  close(con);
  con <- NULL;

  # Rename temporary file
  pathname <- popTemporaryFile(pathnameT, verbose=verbose);

  # Object to be returned
  res <- newInstance(static, pathname);

  # Write footer
  writeFooter(res, footer);

  # Return
  res;
}, static=TRUE, protected=TRUE)



setMethodS3("getColClasses", "AromaTabularBinaryFile", function(this, ...) {
  hdr <- readHeader(this)$dataHeader;
  hdr$types;
})

setMethodS3("getBytesPerColumn", "AromaTabularBinaryFile", function(this, ...) {
  hdr <- readHeader(this)$dataHeader;
  hdr$sizes;
})


setMethodS3("nbrOfRows", "AromaTabularBinaryFile", function(this, ...) {
  hdr <- readHeader(this)$dataHeader;
  hdr$nbrOfRows;
})

setMethodS3("nbrOfColumns", "AromaTabularBinaryFile", function(this, ...) {
  hdr <- readHeader(this)$dataHeader;
  hdr$nbrOfColumns;
})



setMethodS3("[", "AromaTabularBinaryFile", function(this, i=NULL, j=NULL, drop=FALSE) {
  # Read data
  data <- readDataFrame(this, rows=i, columns=j);

  # Drop dimensions?
  if (drop) {
    if (ncol(data) == 1) {
      data <- data[,1];
    } else if (nrow(data) == 1) {
      data <- data[1,];
    }
  }

  data;
})


setMethodS3("[[", "AromaTabularBinaryFile", function(this, i) {
  if (!is.numeric(i))
    throw("Argument 'i' must be numeric: ", i);

  if (length(i) != 1)
    throw("Argument 'i' must be a single value: ", length(i));

  readDataFrame(this, columns=i)[[1]];
}, protected=TRUE)



setMethodS3("[<-", "AromaTabularBinaryFile", function(this, i=NULL, j=NULL, value) {
  updateData(this, rows=i, columns=j, values=value);
  invisible(this);
})


setMethodS3("subset", "AromaTabularBinaryFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  data <- readDataFrame(this);
  subset(data, ...);
})


setMethodS3("summary", "AromaTabularBinaryFile", function(object, ...) {
  # To please R CMD check
  this <- object;

  nbrOfColumns <- nbrOfColumns(this);

  # Get the summaries (as matrices; less work for us, more for R)
  res <- lapply(seq_len(nbrOfColumns), FUN=function(cc) {
    s <- summary(this[,cc,drop=FALSE], ...);
  })

  if (nbrOfColumns == 1) {
    return(res[[1]]);
  }

  # Get the summaries (as matrices; less work for us, more for R)
  res <- lapply(res, FUN=function(s) {
    dimnames(s) <- NULL;
    s <- strsplit(s, split=":");
    names <- lapply(s, FUN=function(str) str[1]);
    values <- lapply(s, FUN=function(str) str[2]);
    names(values) <- names;
    values;
  })

  names <- lapply(res, FUN=function(s) names(s));
  unames <- unique(unlist(names, use.names=FALSE));
  emptyName <- paste(rep(" ", nchar(unames[1])+1), collapse="");

  for (kk in seq_along(res)) {
    s <- res[[kk]];
    emptyStr <- paste(rep(" ", nchar(s[[1]])), collapse="");
    thisNames <- names[[kk]];
    idx <- match(unames, thisNames);
    s <- s[idx];
    nok <- which(is.na(idx));
    s[nok] <- emptyStr;
    thisNames <- paste(thisNames, ":", sep="");
    thisNames[nok] <- emptyName;
    s <- paste(thisNames, s, sep="");
    res[[kk]] <- s;
  }

  res <- matrix(unlist(res, use.names=FALSE), ncol=nbrOfColumns);
  rownames(res) <- rep("", nrow(res));
  colnames(res) <- getColumnNames(this);
  class(res) <- "table";

  res;
})


setMethodS3("colApply", "AromaTabularBinaryFile", function(X, FUN, ...) {
  # To please R CMD check
  this <- X;

  # Argument 'FUN':
  FUN <- match.fun(FUN);

  nbrOfColumns <- nbrOfColumns(this);
  res <- lapply(seq_len(nbrOfColumns), FUN=function(cc) {
    FUN(this[[cc]], ...);
  });

  res;
}, protected=TRUE)


setMethodS3("colStats", "AromaTabularBinaryFile", function(this, FUN, ...) {
  res <- colApply(this, FUN=FUN, ...);
  res <- unlist(res, use.names=FALSE);
  res;
}, protected=TRUE)


setMethodS3("colSums", "AromaTabularBinaryFile", function(x, ...) {
  colStats(x, FUN=sum, ...);
})

setMethodS3("colMeans", "AromaTabularBinaryFile", function(x, ...) {
  colStats(x, FUN=mean, ...);
})

setMethodS3("importFrom", "AromaTabularBinaryFile", function(this, srcFile, ...) {
  methodNames <- sprintf("importFrom%s", class(srcFile));

  # Search for importFrom<ClassName>() methods
  keep <- sapply(methodNames, FUN=exists, mode="function");
  methodNames <- methodNames[keep];

  # Failure?
  if (length(methodNames) == 0) {
    throw("Cannot import from ", class(srcFile)[1], ". Failed to locate importFrom<ClassName>() method.");
  }

  # Use the first method
  methodName <- methodNames[1];
  fcn <- get(methodName, mode="function");

  # Import data
  fcn(this, srcFile=srcFile, ...);
})


############################################################################
# HISTORY:
# 2015-05-03
# o Relaxed sanity checks; now it's possible to allocate Aroma tabular
#   binary files with 200e6 rows (was 100e6 rows).
# 2012-11-29
# o Renamed lapply() for AromaTabularBinaryFile to colApply().
# 2011-09-24
# o readHeader(), readRawFooter() and writeRawFooter() of
#   AromaTabularBinaryFile would try to read non-signed 4-byte integers,
#   which is not supported and would be instead read as signed integers.
#   From R v2.13.1 this would generated warnings.
# 2011-02-01
# o ROBUSTNESS: Using argument 'nchars' (not 'nchar') in readChar() calls.
# 2010-06-02
# o BUG FIX: updateDataColumn() of AromaTabularBinaryFile would
#   censor *signed integers* incorrectly; it should censor at/to
#   [-(n+1),n], but did it at [-n,(n+1)] ("two's complement").
#   This caused it to write too large values as n+1, which then
#   would be read as -(n+1), e.g. writing 130 would be censored
#   to 128 (should be 127), which then would be read as -128.
#   Added more detailed information on how many values were censored.
#   Thanks Robert Ivanek for report on this.
# 2010-01-06
# o Added argument 'defaults' to allocate() of AromaTabularBinaryFile.
# 2009-05-18
# o All methods for modifying an existing file, are now calling
#   Arguments$getWritablePathname() to assert that it exist and that
#   the permissions allow it to be modified.
# 2009-05-16
# o ROBUSTNESS: Now allocate() for AromaTabularBinaryFile first allocates
#   a temporary file which is then renamed.  This makes the file allocation
#   more atomic, and therefore less error prone to interrupts etc.
# o ROBUSTNESS: Now all "write" methods that updates a file asserts that
#   the file exists with isFile(pathname) before open it with
#   file(pathname, open="r+b").
# 2009-05-03
# o Now readDataFrame() of AromaTabularBinaryFile accepts rows=integer(0).
# 2008-12-03
# o CLEAN UP: readDataFrame() of AromaTabularBinaryFile would forget to
#   close the connection if verbose output was activated.  When R later
#   would close such connections, a warning would be generated.
# 2008-07-21
# o Now updateDataColumn() coerce values to doubles before censoring them
#   for raw and integer columns.
# 2008-07-10
# o Added argument 'drop=FALSE' to readDataFrame().
# o SPEED UP: Removed gc() in inner loop of readDataFrame().
# o Now readColumns() calls readDataFrame() and not vice versa.
# 2008-07-09
# o Added a general importFrom() that calls importFrom<ClassName>().
# o Added support for 'raw' columns.
# 2008-05-25
# o BUG FIX: In several methods, when using verbose the on.exit() to close
#   an opened connection was overwritten.
# 2008-05-21
# o TYPO: File-format error message for incorrect 'magic' sequence was
#   reporting the incorrect true sequence.
# o Now allocate() takes argument 'footer' as well.
# 2008-05-16
# o Added readColumns().
# o Inherits from GenericTabularFile.
# 2008-05-11
# o Added extractMatrix() with returns a matrix with one column.  This
#   will be used by the corresponding method for the file set.
# 2008-04-14
# o Renamed readData() to readDataFrame() for AromaTabularBinaryFile.
# 2008-02-13
# o Added and updated Rdoc comments.
# o readFooter() and writeFooter() now handles nested XML/list structures.
# 2008-01-19
# o Now writeFooter() writes a named list to the binary file footer as a
#   XML string.  readFooter() reads its back and returns a named list.
# o Added writeRawFooter().
# o BUG FIX: readRawFooter() did not locate the footer correctly.
# 2007-09-14
# o Renamed static method create() to allocate().
# 2007-09-13
# o Added a file footer to the file format, which comes after the data
#   section.  This way we can append any amount of meta data to the file
#   after it has been created/allocated.
# o Now integers out of range are censored to the limits with a warning.
# 2007-09-11
# o Added colStats(), colSums(), colMeans(), and colMedians().
# o Added dim(), nrow(), ncol() and a memory efficient lapply().
# o Created from AromaGenomeInformationFile.R.  This file type is more
#   general.  There will be an option for a file header too, or rather a
#   file tail, because that is easier to expand.
############################################################################
