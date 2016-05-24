setConstructorS3("AromaUnitChromosomeTabularBinaryFile", function(...) {
  this <- extend(AromaUnitTabularBinaryFile(...), "AromaUnitChromosomeTabularBinaryFile",
    "cached:.memoryCache" = list(),
    .chromosomes=NULL
  );

  # Parse attributes (all subclasses must call this in the constructor).
  setAttributesByTags(this)

  this;
}, abstract=TRUE)



setMethodS3("getGenomeVersion", "AromaUnitChromosomeTabularBinaryFile", function(this, ...) {
  tags <- getTags(this, ...);
  tags <- grep("^hg", tags, value=TRUE);
  tags;
}, protected=TRUE)



setMethodS3("getFilenameExtension", "AromaUnitChromosomeTabularBinaryFile", abstract=TRUE);


setMethodS3("getDefaultColumnNames", "AromaUnitChromosomeTabularBinaryFile", abstract=TRUE);


setMethodS3("indexOfColumn", "AromaUnitChromosomeTabularBinaryFile", function(this, name, ...) {
  cc <- which(getColumnNames(this) == name);
  cc <- Arguments$getIndex(cc);
  cc;
}, protected=TRUE)


setMethodS3("getChromosomes", "AromaUnitChromosomeTabularBinaryFile", function(this, force=FALSE, .chromosomes=NULL, ...) {
  chromosomes <- this$.chromosomes;
  if (force || is.null(chromosomes)) {
    chromosomes <- .chromosomes;
    if (is.null(chromosomes)) {
      cc <- indexOfColumn(this, "chromosome");
      # Sanity check
      if (length(cc) == 0) {
        throw(sprintf("Failed to infer set of chromosomes. There is no column 'chromosome' in this %s file: %s", class(this)[1], getPathname(this)));
      }
      chromosomes <- this[,cc,drop=TRUE];
    }
    chromosomes <- unique(chromosomes);
    chromosomes <- chromosomes[!is.na(chromosomes)];
    chromosomes <- sort(chromosomes);
    this$.chromosomes <- chromosomes;
  }

  chromosomes;
})


setMethodS3("readDataFrame", "AromaUnitChromosomeTabularBinaryFile", function(this, rows=units, ..., units=NULL, verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  data <- NextMethod("readDataFrame", rows=rows, verbose=less(verbose));

  if (nrow(data) > 0) {
    verbose && enter(verbose, "Converting zeros to NAs");
    # Interpret zeros as NAs
    if (ncol(data) > 0) {
      for (cc in seq_len(ncol(data))) {
        nas <- (!is.na(data[,cc]) & (data[,cc] == 0));
        data[nas,cc] <- NA;
      }
    }
    verbose && exit(verbose);
  }

  data;
})



setMethodS3("getUnitsOnChromosomes", "AromaUnitChromosomeTabularBinaryFile", function(this, chromosomes=getChromosomes(this), ..., unlist=TRUE, useNames=!unlist) {
  # Argument 'chromosomes':
  chromosomes <- Arguments$getIndices(chromosomes);

  # Argument 'unlist':
  unlist <- Arguments$getLogical(unlist);

  # Argument 'useNames':
  useNames <- Arguments$getLogical(useNames);

  # Stratify by chromosome
  cc <- indexOfColumn(this, "chromosome");
  data <- this[,cc,drop=TRUE];

  # Update known chromosomes, if not already done.
  allChromosomes <- getChromosomes(this, .chromosomes=data);

  res <- vector("list", length(chromosomes));
  for (cc in seq_along(chromosomes)) {
    units <- which(data == chromosomes[cc]);
    res[[cc]] <- units;
  } # for (cc ...)

  if (useNames) {
    names(res) <- sprintf("Chr%02d", chromosomes);
  }

  if (unlist) {
    res <- unlist(res, use.names=useNames);
  }

  # CONTRACT/Sanity check
  if (unlist) {
    res <- Arguments$getIndices(res);
  } else {
    # Ignored; to expensive
  }

  res;
}, protected=TRUE)


setMethodS3("getUnitsOnChromosome", "AromaUnitChromosomeTabularBinaryFile", function(this, chromosome, ...) {
  # Argument 'chromosome':
  chromosome <- Arguments$getIndex(chromosome);

  units <- getUnitsOnChromosomes(this, chromosomes=chromosome,
                                 unlist=TRUE, useNames=FALSE);
  units;
}, protected=TRUE)


setMethodS3("extractByChromosome", "AromaUnitChromosomeTabularBinaryFile", function(this, chromosomes=getChromosomes(this), ...) {
  unitsList <- getUnitsOnChromosomes(this, chromosomes=chromosomes, unlist=FALSE);
  data <- readDataFrame(this, ...);
  data <- cbind(unit=seq_len(nrow(data)), data);
  lapply(unitsList, FUN=function(units) {
    data[units,,drop=FALSE];
  });
}, protected=TRUE)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# BEGIN: File I/O
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("allocate", "AromaUnitChromosomeTabularBinaryFile", function(static, ..., platform, chipType, footer=list()) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'platform':
  platform <- Arguments$getCharacter(platform, length=c(1,1));

  # Argument 'chipType':
  chipType <- Arguments$getCharacter(chipType, length=c(1,1));

  # Argument 'footer':
  if (is.null(footer)) {
  } else if (!is.list(footer)) {
    throw("Argument 'footer' must be NULL or a list: ", class(footer)[1]);
  }


  # Create file footer
  footer <- c(
    list(
      createdOn=format(Sys.time(), "%Y%m%d %H:%M:%S", usetz=TRUE),
      platform=platform,
      chipType=chipType
    ),
    footer
  );

  # Allocate file
  NextMethod("allocate", footer=footer);
}, static=TRUE, protected=TRUE)




############################################################################
# HISTORY:
# 2012-11-08
# o BUG FIX: readDataFrame() must accept argument 'rows', so in order
#   to also support alias argument 'units', we now do 'rows=units' and
#   'units=NULL'.
# 2012-10-31
# o Added argument 'units' to readDataFrame() for
#   AromaUnitChromosomeTabularBinaryFile.
# 2011-03-03
# o ROBUSTNESS: Added a return contract/sanity check asserting that
#   getUnitsOnChromosomes() truly returns valid 'unit' indices.
# 2010-01-25
# o ROBUSTNESS: Added a sanity check getChromosomes() for class
#   AromaUnitChromosomeTabularBinaryFile validating that the file has a
#   'chromosome' column.
# 2009-09-07
# o Now getUnitsOnChromosomes() returns a vector by default (unlist=TRUE).
# o By default, getUnitsOnChromosomes() now returns names if the return
#   structure is a list, otherwise not.
# 2009-05-08
# o Added allocate() to AromaUnitChromosomeTabularBinaryFile.  This will
#   enforce all subclasses to specify a platform and chip type.
# o Added extractByChromosome().
# o Added indexOfColumn().
# o Extracted AromaUnitChromosomeTabularBinaryFile from AromaUgpFile, where
#   the latter now inherits from the former.
############################################################################
