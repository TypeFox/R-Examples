setConstructorS3("UflSnpInformation", function(..., .ufl=NULL, .verify=TRUE) {
  this <- extend(SnpInformation(...), "UflSnpInformation",
    .ufl = .ufl
  )
  if (.verify && isFile(this)) verify(this)
  this
})

setMethodS3("getAromaUflFile", "UflSnpInformation", function(this, ..., force=FALSE) {
  ufl <- this$.ufl;

  if (force || is.null(ufl)) {
    ufl <- AromaUflFile(getPathname(this), ...);
    this$.ufl <- ufl;
  }

  # Sanity check
  ufl <- Arguments$getInstanceOf(ufl, "AromaUflFile");

  ufl;
}, protected=TRUE);


setMethodS3("getChipType", "UflSnpInformation", function(this, ...) {
  ufl <- getAromaUflFile(this, ...);
  chipType <- getChipType(ufl, ...);
  chipType;
})


setMethodS3("findByChipType", "UflSnpInformation", function(static, ...) {
  AromaUflFile$findByChipType(...);
}, static=TRUE, protected=TRUE)


setMethodS3("byChipType", "UflSnpInformation", function(static, chipType, tags=NULL, nbrOfUnits=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'chipType':
  chipType <- Arguments$getCharacter(chipType, length=c(1,1));

  # Argument 'nbrOfUnits':
  if (!is.null(nbrOfUnits)) {
    nbrOfUnits <- Arguments$getInteger(nbrOfUnits, range=c(0,Inf));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  ufl <- AromaUflFile$byChipType(chipType, tags=tags, nbrOfUnits=nbrOfUnits, ...);
  pathname <- getPathname(ufl);

  verbose && enter(verbose, "Instantiating ", class(static)[1]);
  verbose && cat(verbose, "Pathname: ", pathname);
  verbose && cat(verbose, "Arguments:");
  verbose && str(verbose, list(...));

  res <- newInstance(static, filename=pathname, path=NULL, .ufl=ufl,
                                                  .verify=FALSE, ...);
  verbose && print(verbose, res);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validation?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(nbrOfUnits)) {
    if (nbrOfUnits(res) != nbrOfUnits) {
      throw("The number of units in the loaded ", class(static)[1], " does not match the expected number: ", nbrOfUnits(res), " != ", nbrOfUnits);
    }
  }

  verbose && exit(verbose);

  res;
}, static=TRUE)


setMethodS3("verify", "UflSnpInformation", function(this, ...) {
  tryCatch({
    df <- readDataFrame(this, nrow=10);
  }, error = function(ex) {
    throw("File format error of the UFL SNP information file (",
                                 ex$message, "): ", getPathname(this));
  })
  invisible(TRUE);
}, private=TRUE)


setMethodS3("readDataFrame", "UflSnpInformation", function(this, nrow=NULL, ..., verbose=FALSE) {
  verbose && enter(verbose, "Reading data from UFL file");

  ufl <- getAromaUflFile(this);
  verbose && print(verbose, ufl, level=-20);

  if (is.null(nrow)) {
    verbose && cat(verbose, "Reading all ", nbrOfUnits(ufl), " units");
    res <- ufl[,,drop=FALSE];
  } else {
    units <- 1:nrow;
    verbose && cat(verbose, "Reading ", length(units), " units");
    res <- ufl[units,,drop=FALSE];
  }

  colnames(res) <- c("fragmentLength");
  verbose && exit(verbose);

  res;
})

setMethodS3("getDataColumns", "UflSnpInformation", function(this, ...) {
  ufl <- getAromaUflFile(this);
  names <- getColumnNames(ufl);
  names <- gsub("^length", "fragmentLength", names);
  names;
}, private=TRUE)

setMethodS3("getFields", "UflSnpInformation", function(this, ...) {
  getDataColumns(this, ...);
}, protected=TRUE)

setMethodS3("nbrOfUnits", "UflSnpInformation", function(this, ...) {
  ufl <- getAromaUflFile(this);
  nbrOfUnits(ufl);
})


setMethodS3("isCompatibleWithCdf", "UflSnpInformation", function(this, cdf, ...) {
  # Argument 'cdf':
  cdf <- Arguments$getInstanceOf(cdf, "AffymetrixCdfFile");

  res <- FALSE;

  if (nbrOfUnits(this) != nbrOfUnits(cdf)) {
    attr(res, "reason") <- sprintf("The number of units of the %s and the %s does not match: %s != %s", class(this)[1], class(cdf)[1], nbrOfUnits(this), nbrOfUnits(cdf));
    return(res);
  }

  TRUE;
})


setMethodS3("getData", "UflSnpInformation", function(this, units=NULL, fields=getDataColumns(this), orderBy=NULL, ..., force=FALSE, verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  data <- this$.data;
  if (is.null(data) || force) {
    verbose && enter(verbose, "Retrieving SNP information from file");

    # Now read the genome information data
    ufl <- getAromaUflFile(this);
    cc <- match(fields, getDataColumns(this));
    missing <- fields[is.na(cc)];
    if (length(missing)) {
      throw("Unknown fields: ", paste(missing, collapse=", "));
    }

    verbose && enter(verbose, "Reading SNP information data");
    data <- ufl[,,drop=FALSE];
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


setMethodS3("nbrOfEnzymes", "UflSnpInformation", function(this, ...) {
  cols <- getDataColumns(this);
  length(cols);
})

setMethodS3("getFragmentLengths", "UflSnpInformation", function(this, enzymes=seq_len(nbrOfEnzymes(this)), ...) {
  data <- getData(this, ..., fields=getDataColumns(this)[enzymes]);
  fl <- data[,enzymes,drop=FALSE];
  fl <- as.matrix(fl);
  dim <- dim(fl);
  fl <- as.integer(fl);
  dim(fl) <- dim;
  fl;
})

setMethodS3("getFragmentStarts", "UflSnpInformation", function(this, ...) {
  throw("Not supported.");
})


setMethodS3("getFragmentStops", "UflSnpInformation", function(this, ...) {
  throw("Not supported.");
})


############################################################################
# HISTORY:
# 2009-02-10
# o Added optional validation of number of units to byChipType().
# o Static method byChipType() was not declared static.
# 2008-07-23
# o Now isCompatibleWithCdf() adds attribute 'reason' to FALSE explaining
#   why the object is not compatible.
# 2008-01-20
# o Made argument 'chipType' and 'tags' explicit for fromChipType().
# 2007-12-10
# o Added getChipType() to UflSnpInformation.
# 2007-11-19
# o Added nbrOfEnzymes().
# 2007-09-16
# o BUG FIX: getFragmentLengths() of UflSnpInformation would thrown an error
#   reporting "Unknown fields: fragmentLength".  Now getDataColumns()
#   returns the correct names.
# 2007-09-11
# o Created from DChipSnpInformation.R.
############################################################################
