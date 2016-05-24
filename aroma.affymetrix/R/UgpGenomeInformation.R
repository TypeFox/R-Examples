###########################################################################/**
# @RdocClass UgpGenomeInformation
#
# @title "The UgpGenomeInformation class"
#
# \description{
#  @classhierarchy
#
#  This class represents Aroma UGP genome information files.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "GenomeInformation".}
#   \item{.ugp}{For internal use only.}
#   \item{.verify}{For internal use only.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("UgpGenomeInformation", function(..., .ugp=NULL, .verify=TRUE) {
  this <- extend(GenomeInformation(..., .verify=FALSE), "UgpGenomeInformation",
    "cached:.ugp" = .ugp
  )
  if (.verify && isFile(this)) verify(this)
  this
})


setMethodS3("getAromaUgpFile", "UgpGenomeInformation", function(this, ..., force=FALSE) {
  ugp <- this$.ugp;
  if (force || is.null(ugp)) {
    ugp <- AromaUgpFile(getPathname(this), ...);
    this$.ugp <- ugp;
  }
  ugp;
}, protected=TRUE);


setMethodS3("getChipType", "UgpGenomeInformation", function(this, ...) {
  ugp <- getAromaUgpFile(this, ...);
  chipType <- getChipType(ugp, ...);
  chipType;
})


setMethodS3("findByChipType", "UgpGenomeInformation", function(static, ...) {
  AromaUgpFile$findByChipType(...);
}, static=TRUE, protected=TRUE)


###########################################################################/**
# @RdocMethod byChipType
#
# @title "Defines a UgpGenomeInformation object by chip type"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{chipType}{A @character string.}
#  \item{tags}{...}
#  \item{version}{An optional @character string specifying the version
#    string, if more than one version is available.}
#  \item{...}{Not used.}
#  \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns an @see "UgpGenomeInformation" object.
#  If no file was not found, an error is thrown.
# }
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("byChipType", "UgpGenomeInformation", function(static, chipType, tags=NULL, ..., nbrOfUnits=NULL, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'chipType':
  chipType <- Arguments$getCharacter(chipType, length=c(1,1));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  ugp <- AromaUgpFile$byChipType(chipType, tags=tags, nbrOfUnits=nbrOfUnits, ...);
  pathname <- getPathname(ugp);

  verbose && enter(verbose, "Instantiating ", class(static)[1]);
  verbose && cat(verbose, "Pathname: ", pathname);
  verbose && cat(verbose, "Arguments:");
  verbose && str(verbose, list(...));

  res <- newInstance(static, filename=pathname, path=NULL, .ugp=ugp, ...);
  verbose && print(verbose, res);

  verbose && exit(verbose);

  res;
}, static=TRUE)



setMethodS3("nbrOfUnits", "UgpGenomeInformation", function(this, ...) {
  ugp <- getAromaUgpFile(this);
  nbrOfUnits(ugp);
})

setMethodS3("isCompatibleWithCdf", "UgpGenomeInformation", function(this, cdf, ...) {
  # Argument 'cdf':
  cdf <- Arguments$getInstanceOf(cdf, "AffymetrixCdfFile");

  res <- FALSE;

  if (nbrOfUnits(this) != nbrOfUnits(cdf)) {
    attr(res, "reason") <- sprintf("The number of units of the %s and the %s does not match: %s != %s", class(this)[1], class(cdf)[1], nbrOfUnits(this), nbrOfUnits(cdf));
    return(res);
  }

  TRUE;
})



setMethodS3("verify", "UgpGenomeInformation", function(this, ...) {
  tryCatch({
    df <- readDataFrame(this, nrow=10);
  }, error = function(ex) {
    throw("File format error of the UGP genome information file (",
                                 ex$message, "): ", getPathname(this));
  })
  invisible(TRUE);
}, private=TRUE)


setMethodS3("readDataFrame", "UgpGenomeInformation", function(this, units=NULL, nrow=NULL, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Reading data from UGP file");
  ugp <- getAromaUgpFile(this);
  verbose && print(verbose, ugp, level=-20);

  # Validate 'units'
  if (is.null(units)) {
    if (is.null(nrow)) {
      units <- 1:nbrOfUnits(ugp);
    } else {
      units <- 1:nrow;
    }
  }
  units <- Arguments$getIndices(units, max=nbrOfUnits(ugp));


  verbose && enter(verbose, "Reading ", length(units), " units");
  res <- ugp[units,,drop=FALSE];
  verbose && exit(verbose);

  colnames(res) <- c("chromosome", "physicalPosition");
  verbose && str(verbose, res);
  verbose && exit(verbose);

  res;
})


setMethodS3("getDataColumns", "UgpGenomeInformation", function(this, ...) {
  c("chromosome", "physicalPosition");
}, private=TRUE)


setMethodS3("getData", "UgpGenomeInformation", function(this, units=NULL, fields=getDataColumns(this), orderBy=NULL, ..., force=FALSE, verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  data <- this$.data;
  if (is.null(data) || force) {
    verbose && enter(verbose, "Retrieving genome information from file");

    # Now read the genome information data
    ugp <- getAromaUgpFile(this);
    cc <- match(fields, getDataColumns(this));
    missing <- fields[is.na(cc)];
    if (length(missing)) {
      throw("Unknown fields: ", paste(missing, collapse=", "));
    }

    verbose && enter(verbose, "Reading genome information data");
    data <- ugp[,,drop=FALSE];
    colnames(data) <- getDataColumns(this);
    verbose && str(verbose, data);
    verbose && exit(verbose);

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
    data <- data[units,,drop=FALSE];
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


setMethodS3("getChromosomes", "UgpGenomeInformation", function(this, force=FALSE, ...) {
  ugp <- getAromaUgpFile(this);
  getChromosomes(ugp, force=force, ...);
})

setMethodS3("getUnitsOnChromosome", "UgpGenomeInformation", function(this, ...) {
  ugp <- getAromaUgpFile(this);
  getUnitsAt(ugp, ...);
})

############################################################################
# HISTORY:
# 2009-02-16
# o Removed argument 'validate' from byChipType().
# 2008-12-29
# o Added argument 'units' to readDataFrame().
# 2008-07-23
# o Now isCompatibleWithCdf() adds attribute 'reason' to FALSE explaining
#   why the object is not compatible.
# 2008-05-20
# o Added argument 'validate' to byChipType().
# 2008-01-19
# o ROBUSTNESS: Added argument 'tags' to fromChipType() for class
#   UgpGenomeInformation.
# 2007-12-10
# o Added getChipType() to UgpGenomeInformation.
# 2007-12-09
# o Added nbrOfUnits().
# 2007-11-20
# o Added clearCache() to clear cached UGP file.
# 2007-09-11
# o Created.
############################################################################
