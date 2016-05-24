# @author "HB"
setConstructorS3("AffymetrixCsvGenomeInformation", function(...) {
  this <- extend(GenomeInformation(...), "AffymetrixCsvGenomeInformation")
  if (isFile(this)) verify(this)
  this
})


setMethodS3("getDefaultExtension", "AffymetrixCsvGenomeInformation", function(static, ...) {
  "annot.csv";
}, static=TRUE, protected=TRUE)


setMethodS3("findByChipType", "AffymetrixCsvGenomeInformation", function(static, chipType, version=NULL, ...) {
  # Argument 'version':
  if (is.null(version))
    version <- ".*";

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Search in annotationData/chipTypes/<chipType>/
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ext <- getDefaultExtension(static);
  pattern <- sprintf("^%s[.].*[.]%s$", chipType, ext);
  pathname <- findAnnotationDataByChipType(chipType, pattern);

  pathname;
}, static=TRUE, protected=TRUE)


setMethodS3("byChipType", "AffymetrixCsvGenomeInformation", function(static, chipType, version=NULL, ...) {
  # Search for the genome information file
  pathname <- findByChipType(static, chipType, version=version, ...);
  if (is.null(pathname))
    throw("Failed to located Affymetrix CSV annotation file: ", chipType);
  newInstance(static, pathname);
}, static=TRUE)

setMethodS3("verify", "AffymetrixCsvGenomeInformation", function(this, ...) {
  tryCatch({
    df <- readDataFrame(this, nrow=10);
  }, error = function(ex) {
    throw("File format error of the Affymetrix CSV annotation file: ",
                                                  getPathname(this));
  })
  invisible(TRUE);
}, protected=TRUE)


setMethodS3("readDataFrame", "AffymetrixCsvGenomeInformation", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  ext <- getDefaultExtension(this);
  verbose && enter(verbose, sprintf("Reading Affymetrix *.%s file", ext));

  pathname <- getPathname(this);

  verbose && cat(verbose, "Pathname: ", pathname);

  hdr <- scan(pathname, what=character(0), nlines=1, sep=",",
                                           quote="\"", quiet=TRUE);
  nbrOfColumns <- length(hdr);

  verbose && printf(verbose, "Columns [%d]:\n", nbrOfColumns);
  verbose && print(verbose, hdr);

  colClasses <- rep("NULL", nbrOfColumns);
  names(colClasses) <- hdr;

  fields <- c("Probe Set ID", "PROBESET_ID");
  cc <- na.omit(match(fields, names(colClasses)));
  colClasses[cc] <- "character";

  fields <- c("Chromosome", "CHROMOSOME");
  cc <- na.omit(match(fields, names(colClasses)));
  colClasses[cc] <- "character";

  fields <- c("Physical Position", "PROBE_START_POSITION");
  cc <- na.omit(match(fields, names(colClasses)));
  colClasses[cc] <- "character";

  verbose && printf(verbose, "colClasses [%d]:\n", length(colClasses));
  verbose && str(verbose, as.list(colClasses));

  # Make sure we haven't added or removed columns
  stopifnot(length(colClasses) == nbrOfColumns);

  # Read the data table
  df <- read.table(pathname, colClasses=colClasses, header=TRUE, sep=",", quote="\"", fill=TRUE, check.names=FALSE, na.strings=c("---"), ...);

  # Update the column names
  colnames <- colnames(df);

  fields <- c("Probe Set ID", "PROBESET_ID");
  cc <- na.omit(match(fields, colnames));
  colnames[cc] <- fields[1];

  fields <- c("Chromosome", "CHROMOSOME");
  cc <- na.omit(match(fields, colnames));
  colnames[cc] <- fields[1];

  fields <- c("Physical Position", "PROBE_START_POSITION");
  cc <- na.omit(match(fields, colnames));
  colnames[cc] <- fields[1];

  colnames <- toCamelCase(colnames);
  colnames(df) <- colnames;

  # Chromosome
  chr <- df[["chromosome"]];
  chr[chr == "X"] <- 23;
  chr[chr == "Y"] <- 24;
  suppressWarnings({
    chr <- as.integer(chr);
  })
  df[["chromosome"]] <- chr;
  # Not needed anymore
  chr <- NULL;

  # Coerce to integers
  df[["physicalPosition"]] <- as.integer(df[["physicalPosition"]]);

  verbose && exit(verbose);

  df;
})


############################################################################
# HISTORY:
# 2011-11-19
# o Added getDefaultExtension() for AffymetrixCsvGenomeInformation.
# 2008-04-14
# o Renamed readData() to readDataFrame() for AffymetrixCsvGenomeInformation.
# 2007-08-03
# o BUG FIX: Declared fromChipType() static.
# 2007-03-02
# o Created.
############################################################################
