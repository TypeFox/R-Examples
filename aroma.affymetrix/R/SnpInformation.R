###########################################################################/**
# @RdocClass SnpInformation
#
# @title "The SnpInformation class"
#
# \description{
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "R.filesets::GenericDataFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("SnpInformation", function(...) {
  this <- extend(GenericDataFile(...), c("SnpInformation"
                                             , uses("FileCacheKeyInterface")),
    "cached:.data" = NULL
  )
  if (isFile(this)) verify(this)
  this
})

setMethodS3("as.character", "SnpInformation", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character");
  s <- c(s, sprintf("Chip type: %s", getChipType(this)));
  s <- c(s, sprintf("Number of enzymes: %s", nbrOfEnzymes(this)));
  s;
}, protected=TRUE)



###########################################################################/**
# @RdocMethod verify
#
# @title "Verifies the correctness of the underlying file"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns (visibly) @TRUE if the file is valid, otherwise an error is
#   thrown.
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("verify", "SnpInformation", function(this, ...) {
  TRUE;
}, protected=TRUE)



###########################################################################/**
# @RdocMethod getChipType
#
# @title "Gets the chip type of this genome information set"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character string.
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getChipType", "SnpInformation", function(this, ...) {
  # Infer chip type from the first parent directory that has the same name
  # as the chip type of an existing CDF file.
  pathname <- getPathname(this);
  lastPath <- pathname;
  while (TRUE) {
    path <- dirname(lastPath);
    if (path == lastPath)
      break;
    chipType <- basename(path);
    dummy <- AffymetrixCdfFile$findByChipType(chipType);
    if (!is.null(dummy))
      return(chipType)
    lastPath <- path;
  }
  throw("Failed to infer the chip type from the pathname of the SNP information file: ", pathname);
})


###########################################################################/**
# @RdocMethod fromCdf
#
# @title "Static method to define a genome information set from a CDF"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{cdf}{An @see "AffymetrixCdfFile".}
#  \item{...}{Additional argument passed to @seemethod "byChipType".}
# }
#
# \value{
#   Returns a @see "SnpInformation" object.
# }
#
# \seealso{
#   @seemethod "byChipType".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("fromCdf", "SnpInformation", function(static, cdf, ...) {
  byChipType(static, chipType=getChipType(cdf), nbrOfUnits=nbrOfUnits(cdf), ...);
}, static=TRUE, protected=TRUE)




###########################################################################/**
# @RdocMethod byChipType
#
# @title "Static method to define a genome information set by chip type"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @see "SnpInformation" object.
# }
#
# \seealso{
#   @seemethod "fromCdf".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("byChipType", "SnpInformation", abstract=TRUE);


setMethodS3("fromDataSet", "SnpInformation", function(static, dataSet, ...) {
  chipType <- getChipType(dataSet);
  byChipType(static, chipType=chipType, ...);
}, static=TRUE, protected=TRUE)



setMethodS3("isCompatibleWithCdf", "SnpInformation", function(this, cdf, ...) {
  # Argument 'cdf':
  cdf <- Arguments$getInstanceOf(cdf, "AffymetrixCdfFile");

  # By default, be naive and always return FALSE.
  TRUE;
}, protected=TRUE)



###########################################################################/**
# @RdocMethod getData
#
# @title "Gets all or a subset of the genome information data"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{units}{The units for which the data should be returned.}
#  \item{fields}{The fields to be returned.}
#  \item{orderBy}{The fields by which the returned data frame should be
#      ordered.}
#  \item{...}{Named arguments used to select a subset of the units to be
#      returned.  Either a value to be compared to or a @function returning
#      @TRUE or @FALSE.}
# }
#
# \value{
#   Returns a @data.frame, where the row names correspond to unit indices
#   as defined by the CDF.
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getData", "SnpInformation", function(this, units=NULL, fields=c("fragmentLength", "start", "stop"), orderBy=NULL, ..., force=FALSE, verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  data <- this$.data;
  if (force || is.null(data)) {
    verbose && enter(verbose, "Retrieving SNP information from file");

    # Read the unit names from the corresponding CDF file
    verbose && enter(verbose, "Reading unit names from CDF file");
    chipType <- getChipType(this);
    cdfFile <- AffymetrixCdfFile$findByChipType(chipType);
    if (is.null(cdfFile))
      throw("Could not located CDF file: ", chipType);
    targetUnitNames <- .readCdfUnitNames(cdfFile);
    verbose && exit(verbose);

    # Now read the SNP information data
    verbose && enter(verbose, "Reading SNP information data");
    data <- readDataFrame(this, verbose=less(verbose));
    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);
    verbose && exit(verbose);

    verbose && enter(verbose, "Splitting up the fragment length, start & stop details");
    cc <- which("fragmentLengthStartStop" == colnames(data));
    lss <- data[,cc,drop=TRUE];
    nas <- ((lss == "---") | (lss == ""));
    lss[nas] <- "-//-//-";
    lss <- strsplit(lss, split="//");
    len <- sapply(lss, FUN=length);
    if (any(len != 3)) {
      ulen <- unique(len);
      throw("Internal error: Unrecognized length of FLSS: ",
                                           paste(ulen, collapse=", "));
    }
    lss <- unlist(lss, use.names=FALSE);
    lss <- as.integer(lss);
    lss <- matrix(lss, ncol=3, byrow=TRUE);
    rownames(lss) <- 1:nrow(lss);
    colnames(lss) <- c("fragmentLength", "start", "stop");
    data <- cbind(data[,-cc,drop=FALSE], lss);
    # Not needed anymore
    lss <- NULL;
    verbose && exit(verbose);

    verbose && enter(verbose, "Reordering units according to the CDF file");
    idxs <- match(targetUnitNames, data[,1]);
    data <- data[idxs,];
    rownames(data) <- 1:nrow(data);
#    data <- data[idxs,-1];
    verbose && exit(verbose);
    # Store in cache
    this$.data <- data;

    verbose && exit(verbose);
  }

  # Subset by unit?
  if (!is.null(units)) {
    # Map the unit indicies to the row names
    rr <- match(units, rownames(data));
    data <- data[rr,,drop=TRUE];
  }

  # Stratify by field values?
  args <- list(...);
  if (length(args) > 0) {
    for (key in names(args)) {
      # Get the values to be stratified upon.
      values <- data[,key];

      # Get the test (value or function)
      test <- args[[key]];
      test <- na.omit(test);
      if (is.function(test)) {
        keep <- test(values);
      } else {
        keep <- (values == test);
      }
      data <- data[keep,,drop=FALSE];
    }
    # Not needed anymore
    keep <- NULL;
  }

  # Reorder?
  if (!is.null(orderBy)) {
    o <- do.call(order, args=as.list(data[,orderBy]));
    data <- data[o,,drop=FALSE];
    # Not needed anymore
    o <- NULL;
  }

  # Extract a subset of fields?
  if (!is.null(fields))
    data <- data[,fields, drop=FALSE];

  data;
})


setMethodS3("nbrOfUnits", "SnpInformation", function(this, ...) {
  data <- getData(this, fields=1);
  nrow(data);
})

setMethodS3("getFields", "SnpInformation", function(this, ...) {
  data <- getData(this);
  colnames(data);
}, protected=TRUE)



setMethodS3("readDataFrame", "SnpInformation", abstract=TRUE);


setMethodS3("readTableInternal", "SnpInformation", function(this, pathname, colClasses=NULL, ..., include=NULL, exclude=NULL, verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  verbose && enter(verbose, "Reading tabular data from file");
  pathname <- getPathname(this);
  verbose && cat(verbose, "Pathname: ", pathname);

  verbose && cat(verbose, "Argument 'include': ", paste(include, collapse=", "));
  verbose && cat(verbose, "Argument 'exclude': ", paste(exclude, collapse=", "));
  exclude <- setdiff(exclude, include);
  verbose && cat(verbose, "exclude\\include: ", paste(exclude, collapse=", "));
  colClasses[names(colClasses) %in% exclude] <- "NULL";
  toRead <- names(colClasses)[colClasses != "NULL"];
  verbose && cat(verbose, "Columns to be read: ", paste(toRead, collapse=", "));

  df <- readTable(pathname, colClasses=colClasses, header=TRUE, sep="\t", ..., verbose=less(verbose));

  colnames(df) <- toCamelCase(colnames(df));
  verbose && exit(verbose);

  df;
}, private=TRUE)


setMethodS3("nbrOfEnzymes", "SnpInformation", function(this, ...) {
  as.integer(1);
})

setMethodS3("getFragmentLengths", "SnpInformation", function(this, enzymes=seq_len(nbrOfEnzymes(this)), ...) {
  data <- getData(this, ..., fields="fragmentLength");
  fl <- data[,enzymes,drop=FALSE];
  fl <- as.matrix(fl);
  dim <- dim(fl);
  fl <- as.integer(fl);
  dim(fl) <- dim;
  fl;
})


setMethodS3("getFragmentStarts", "SnpInformation", function(this, enzymes=seq_len(nbrOfEnzymes(this)), ...) {
  data <- getData(this, ..., fields="start");
  fl <- data[,enzymes,drop=FALSE];
  fl <- as.matrix(fl);
  dim <- dim(fl);
  fl <- as.integer(fl);
  dim(fl) <- dim;
  fl;
})


setMethodS3("getFragmentStops", "SnpInformation", function(this, enzymes=seq_len(nbrOfEnzymes(this)), ...) {
  data <- getData(this, ..., fields="stop");
  fl <- data[,enzymes,drop=FALSE];
  fl <- as.matrix(fl);
  dim <- dim(fl);
  fl <- as.integer(fl);
  dim(fl) <- dim;
  fl;
})



############################################################################
# HISTORY:
# 2008-04-14
# o Renamed readData() to readDataFrame() for SnpInformation.
# 2007-11-19
# o Now getFragmentLength/Starts/Stops() of SnpInformation return a matrix
#   where each column correspond to an enzyme. This is was added because
#   the new SNP chips have two enzymes. Added nbrOfEnzymes().
# 2007-01-22
# o BUG FIX: getData() did not support the dChip SNP information file for
#   the Mapping10K_Xba142 chip type, because yet again it used a slightly
#   different format compared with the 100K and the 500K files.
# 2006-09-17
# o Created from GenomeInformation.R.
############################################################################
