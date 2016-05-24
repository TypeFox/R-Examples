# @author "HB"
setConstructorS3("AffymetrixTsvFile", function(...) {
  this <- extend(AffymetrixFile(...), "AffymetrixTsvFile",
    "cached:.cdf" = NULL,
    "cached:.data" = NULL
  )
  if (isFile(this)) verify(this)
  this
})

setMethodS3("getChipType", "AffymetrixTsvFile", function(this, ...) {
  getName(this);
})

setMethodS3("getDefaultExtension", "AffymetrixTsvFile", function(static, ...) {
  "tsv";
}, static=TRUE, protected=TRUE)

setMethodS3("getExtensionPattern", "AffymetrixTsvFile", function(static, ...) {
  ext <- getDefaultExtension(static, ...);
  pattern <- sprintf("[.](%s|%s)$", tolower(ext), toupper(ext));
  pattern;
}, static=TRUE, protected=TRUE)


setMethodS3("findByChipType", "AffymetrixTsvFile", function(static, chipType, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Search in annotationData/chipTypes/<chipType>/
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pattern <- sprintf("^%s%s", chipType, getExtensionPattern(static));
  pathname <- findAnnotationDataByChipType(chipType, pattern);

  pathname;
}, static=TRUE, protected=TRUE)


setMethodS3("byChipType", "AffymetrixTsvFile", function(static, chipType, ...) {
  # Search for the genome information file
  pathname <- findByChipType(static, chipType, ...);
  if (is.null(pathname))
    throw("Failed to located Affymetrix TSV file: ", chipType);
  newInstance(static, pathname);
})

setMethodS3("verify", "AffymetrixTsvFile", function(this, ...) {
  tryCatch({
    df <- readDataFrame(this, nrows=10);
  }, error = function(ex) {
    throw("File format error of the Affymetrix TSV file: ",
                                                  getPathname(this));
  })
  invisible(TRUE);
}, private=TRUE)

setMethodS3("getData", "AffymetrixTsvFile", function(this, force=FALSE, ...) {
  data <- this$.data;
  if (force || is.null(data)) {
    data <- readDataFrame(this, ...);
    this$.data <- data;
  }
  data;
})

setMethodS3("readDataFrame", "AffymetrixTsvFile", function(this, ..., verbose=FALSE) {
  pathname <- getPathname(this);

  colClasses <- c(
    "probeset_id"="character",
    "chr"="character",
    "snp_pos"="integer",
    "len"="integer",
    "GC"="double",
    "gc_count"="integer"
  );

  df <- readTable(pathname, colClasses=colClasses, header=TRUE, sep="\t", ...);

  names <- colnames(df);
  names <- gsub("probeset_id", "unit", names);
  names <- gsub("chr", "chromosome", names);
  names <- gsub("GC", "gc", names);
  names <- gsub("snp_pos", "physical position", names);
  names <- gsub("_", " ", names);
  names <- toCamelCase(names);
  colnames(df) <- names;

  # Rescale GC contents to [0,1]
  df[["gc"]] <- df[["gc"]]/100;

  # Remap chromsome X->23, Y->24
  chr <- df[["chromosome"]];
  chr[chr == "X"] <- 23;
  chr[chr == "Y"] <- 24;
  suppressWarnings({
    chr <- as.integer(chr);
  })
  df[["chromosome"]] <- chr;
  # Not needed anymore
  chr <- NULL;

  # Remove duplicated rows (type by Affymetrix?!? /HB 2007-04-02)
  df <- unique(df);

  gc();

  # Convert unit names to unit indices
  unf <- getUnitNamesFile(this);
  unitNames <- getUnitNames(unf);
  units <- match(df[["unit"]], unitNames);
  if (any(is.na(units))) {
    throw("File format error: Identified units that do not exist in the annotation unit names file: ", getChipType(unf));
  }
  df[["unit"]] <- units;

#  rownames(df) <-  units;

  o <- order(units);
  df <- df[o,];

  df;
})


setMethodS3("getField", "AffymetrixTsvFile", function(this, units=NULL, field, ...) {
  if (is.null(units)) {
    unf <- getUnitNamesFile(this);
    units <- seq_len(nbrOfUnits(unf));
  }

  data <- getData(this, ...);
  if (!field %in% colnames(data))
    throw("No such field: ", field);

  idxs <- match(units, data$unit);
  data[[field]][idxs];
})

setMethodS3("getPosition", "AffymetrixTsvFile", function(this, ...) {
  getField(this, field="physicalPosition", ...);
})

setMethodS3("getFragmentLengths", "AffymetrixTsvFile", function(this, ...) {
  getField(this, field="len", ...);
})

setMethodS3("getGc", "AffymetrixTsvFile", function(this, ...) {
  getField(this, field="gc", ...);
})



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# BEGIN: Methods that requires an Affymetrix CDF file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("getCdf", "AffymetrixTsvFile", function(this, ...) {
  cdf <- this$.cdf;
  if (is.null(cdf)) {
    chipType <- getChipType(this);
    cdf <- AffymetrixCdfFile$byChipType(chipType);
    this$.cdf <- cdf;
  }
  cdf;
})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# END: Methods that requires an Affymetrix CDF file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


############################################################################
# HISTORY:
# 2011-11-19
# o Added getDefaultExtension() to AffymetrixCsvFile.
# 2008-05-18
# o Now readDataFrame() and getField() of AffymetrixTsvFile utilize the
#   UnitNamesFile interface rather than the platform-specific
#   AffymetrixCdfFile.  This is done in order minimize dependencies for
#   certain file formats, i.e. not all chip types comes with a CDF.
# 2008-04-14
# o Renamed readData() to readDataFrame() for AffymetrixTsvFile.
# 2007-03-02
# o Created.
############################################################################
