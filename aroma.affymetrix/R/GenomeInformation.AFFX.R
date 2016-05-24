###########################################################################/**
# @set "class=GenomeInformation"
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
setMethodS3("getChipType", "GenomeInformation", function(this, ...) {
  # Infer chip type from the first parent directory that has the same name
  # as the chip type of an existing annotation unit names file.
  pathname <- getPathname(this);
  lastPath <- pathname;
  while (TRUE) {
    path <- dirname(lastPath);
    if (path == lastPath)
      break;
    chipType <- basename(path);

    # Try to find annotation data file for this chip type
    dummy <- AffymetrixCdfFile$findByChipType(chipType);
    if (!is.null(dummy))
      return(chipType);

    lastPath <- path;
  }

  throw("Failed to infer the chip type from the pathname of the genome information file: ", pathname);
})


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
#   as defined by the annotation unit names file.
# }
#
# \seealso{
#   @seemethod "getUnitIndices".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getData", "GenomeInformation", function(this, units=NULL, fields=c("chromosome", "physicalPosition"), orderBy=NULL, ..., force=FALSE, verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  data <- this$.data;
  if (is.null(data) || force) {
    verbose && enter(verbose, "Retrieving genome information from file");

    # Read the unit names from the corresponding annotation file
    verbose && enter(verbose, "Reading unit names from annotation file");
    chipType <- getChipType(this);
#    unf <- UnitNamesFile$byChipType(chipType);
    cdf <- AffymetrixCdfFile$byChipType(chipType);
    unf <- cdf;
    targetUnitNames <- getUnitNames(unf);
    verbose && exit(verbose);

    # Now read the genome information data
    verbose && enter(verbose, "Reading genome information data");
    data <- readDataFrame(this, verbose=less(verbose));
    verbose && str(verbose, data);
    verbose && exit(verbose);

    verbose && enter(verbose, "Reordering units according to the unit names file");
    idxs <- match(targetUnitNames, data[,1]);
    data <- data[idxs,,drop=FALSE];
#    data <- data[,-1];
    rownames(data) <- 1:nrow(data);
    verbose && exit(verbose);

    verbose && enter(verbose, "Optimizing default return order");
    # Default ordering
    args <- as.list(data[,fields,drop=FALSE]);
    o <- do.call(order, args=args);
    data <- data[o,,drop=FALSE];
    # Not needed anymore
    o <- NULL;
    verbose && str(verbose, data);
    verbose && exit(verbose);

    if ("chromosome" %in% fields) {
      verbose && enter(verbose, "Replacing 'X', 'Y' & 'M' with 23, 24 & 25");
      chr <- data[,"chromosome"];
      chr[chr == "X"] <- 23;
      chr[chr == "Y"] <- 24;
      chr[chr == "M"] <- 25;  # Mitochondrial
      suppressWarnings({
        chr <- as.integer(chr);
      })
      data[,"chromosome"] <- chr;
      # Not needed anymore
      chr <- NULL;
      verbose && str(verbose, data);
      verbose && exit(verbose);
    }

    # Store in cache
    this$.data <- data;

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);

    verbose && exit(verbose);
  }

  # Subset by unit?
  if (!is.null(units)) {
    # Map the unit indicies to the row names
    rr <- match(units, rownames(data));
    data <- data[rr,,drop=FALSE];
  }

  # Stratify by field values?
  args <- list(...);
  if (length(args) > 0) {
    for (key in names(args)) {
      # Get the values to be stratified upon.
      values <- data[,key,drop=FALSE];

      # Get the test (value or function)
      test <- args[[key]];
      test <- na.omit(test);
      if (is.function(test)) {
        keep <- test(values);
      } else {
        keep <- (values == test);
        keep <- (keep & !is.na(keep));
      }
      data <- data[keep,,drop=FALSE];
    }
    # Not needed anymore
    keep <- NULL;
  }

  # Reorder?
  if (!is.null(orderBy)) {
    o <- do.call(order, args=as.list(data[,orderBy,drop=FALSE]));
    data <- data[o,,drop=FALSE];
    # Not needed anymore
    o <- NULL;
  }

  # Extract a subset of fields?
  if (!is.null(fields))
    data <- data[,fields, drop=FALSE];

  data;
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
#   Returns a @see "GenomeInformation" object.
# }
#
# \seealso{
#   @seemethod "byChipType".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("fromCdf", "GenomeInformation", function(static, cdf, ...) {
  byChipType(static, chipType=getChipType(cdf), nbrOfUnits=nbrOfUnits(cdf), ...);
}, static=TRUE, protected=TRUE)



setMethodS3("isCompatibleWithCdf", "GenomeInformation", function(this, cdf, ...) {
  # Argument 'cdf':
  cdf <- Arguments$getInstanceOf(cdf, "AffymetrixCdfFile");

  # By default, be naive and always return...
  TRUE;
}, protected=TRUE)



############################################################################
# HISTORY:
# 2008-05-17
# o Created GenomeInformation.AFFX.R to separate methods specific to
#   Affymetrix from platform generic methods.
# 2008-04-14
# o Renamed readData() to readDataFrame() for GenomeInformation.
# 2007-12-01
# 2007-11-25
# o Now getUnitsOnChromosome() of GenomeInformation can indentify units from
#   multiple chromosomes.
# 2007-10-30
# o Now argument 'chromosome' for getUnitsOnChromosome() needs to be
#   specified explictly. Before its default was '23'.
# 2007-09-10
# o Now readData() in GenomeInformation is no longer abstract, but tries
#   naively to read data using readTableInternal() as is.  That will make
#   it slightly easier to add new genome information files.
# 2007-08-12
# o BUG GIX: clearCache() would not clear the genome stats cache.
# o BUG FIX: Subsetting with getData(..., chromosome=...) would return NAs
#   for units with missing information.
# 2007-06-11
# o BUG FIX: Used non-existing 'gi' instead of 'this' in plotDensity() of
#   GenomeInformation.
# 2007-03-15
# o Updated GenomeInformation to return chromosomes as indices and never
#   with 'X' and 'Y' regardless of source.  This is part of a moving the
#   package to handle chromosomes as indices so that it will work with
#   any genome.
# 2007-02-28
# o Added argument 'region' to getUnitsOnChromosome().
# 2007-01-25
# o BUG FIX: Added 'fill=TRUE' to readTableInternal().
# 2007-01-06
# o Renamed getFields() to getDataColumns().
# 2006-11-29
# o Added getUnitsOnChromosome().
# 2006-09-16
# o Added plotDensity(), getChromosomes(), getChromosomeStats() etc.
# o Improved getData().  Updated getUnitIndices() and getPositions().
# 2006-09-15
# o Added Rdoc comments.
# o Created from DChip.R in old(!) aroma.snp.
# 2005-11-15
# o Added support for skipping header in readSampleInformationFile().
# 2005-10-31
# o Created.
############################################################################

