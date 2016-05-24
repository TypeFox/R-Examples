setMethodS3("readDataFrame", "AffymetrixCdfFile", function(this, units=NULL, fields="*", ..., force=FALSE, verbose=FALSE) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  readCdfDataFrame <- affxparser::readCdfDataFrame


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'units':
  nbrOfUnits <- nbrOfUnits(this);
  if (is.null(units)) {
  } else {
    # Validate unit indices
    units <- Arguments$getIndices(units, max=nbrOfUnits);
    nbrOfUnits <- length(units);
  }

  # Argument 'fields':
  if (!is.null(fields)) {
    fields <- Arguments$getCharacters(fields);
    fields <- unlist(strsplit(fields, split=",", fixed=TRUE));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Reading CDF as a data frame");

  pathname <- getPathname(this);
  verbose && cat(verbose, "Pathname: ", pathname);

  knownVirtualFields <- c("isPm");
  virtualFields <- intersect(fields, knownVirtualFields);
  verbose && cat(verbose, "Virtual fields: ", paste(virtualFields, collapse=", "));

  verbose2 <- as.integer(isVisible(verbose, -10));

##  cdfFields <- setdiff(fields, virtualFields);
##  if (cdfFields == "*") {
##    cdfFields <- NULL;
##  }
##  if ("isPm" %in% virtualFields) {
##    cdfFields <- c(cdfFields, "pbase", "tbase");
##  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Checking for cached results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Checking cached results");
  args <- list(units=units, ...);
  keep <- intersect(names(args), names(formals(readCdfDataFrame)));
  keep <- setdiff(keep, "verbose");
  args <- args[keep];
  key <- list(method="readDataFrame", class=class(this)[1],
              chipType=getChipType(this, fullname=TRUE));
  if (getOption(aromaSettings, "devel/useCacheKeyInterface", FALSE)) {
    key <- getCacheKey(this, method="readDataFrame", chipType=getChipType(this, fullname=TRUE));
  }
  key <- c(key, args);
  verbose && str(verbose, key);
  dirs <- c("aroma.affymetrix", getChipType(this, fullname=TRUE));
  if (force) {
    res <- NULL;
  } else {
    res <- loadCache(key, dirs=dirs);
  }
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Reading data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(res)) {
    verbose && enter(verbose, "Reading data from CDF file");
    verbose && cat(verbose, "Pathname: ", pathname);
    verbose && cat(verbose, "Units:");
    verbose && str(verbose, units);

    t0 <- processTime();
    res <- doCall("readCdfDataFrame", args=list(filename=pathname, units=units, ..., verbose=verbose2));
    t1 <- processTime();

    verbose && cat(verbose, "Data read:");
    verbose && str(verbose, res);

    if (verbose) {
      dt <- (t1-t0)[3];
      cat(verbose, "Total reading/processing time:");
      printf(verbose, "Total time: %.0f secs = %.2f mins = %.2f hours\n",
                                                        dt, dt/60, dt/3600);
      printf(verbose, "Time/unit: %.2f ms = %.2f secs\n",
                                         1000*dt/nbrOfUnits, dt/nbrOfUnits);
    }
    verbose && exit(verbose);

    # Make nucleotide bases in upper case.
    for (field in c("pbase", "tbase")) {
      res[[field]] <- toupper(res[[field]]);
    }

    # Save to file cache
    saveCache(res, key=key, dirs=dirs);
  } else {
    verbose && cat(verbose, "Found results cached on file");
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Generate virtual fields
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdfFields <- names(res);

  if (length(virtualFields) > 0) {
    verbose && enter(verbose, "Adding virtual fields");
    if ("isPm" %in% virtualFields) {
      res[["isPm"]] <- with(res,
        (tbase == "A" & pbase == "T") |
        (tbase == "T" & pbase == "A") |
        (tbase == "C" & pbase == "G") |
        (tbase == "G" & pbase == "C")
      );
    }
    verbose && exit(verbose);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract fields of interest
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(fields)) {
    if (any(fields == "*")) {
      verbose && enter(verbose, "Expanding asterisk fields");
      verbose && cat(verbose, "Argument 'fields': ", paste(fields, collapse=", "));
      asteriskFields <- paste(cdfFields, collapse=",");
      fields[(fields == "*")] <- asteriskFields;
      fields <- unlist(strsplit(fields, split=",", fixed=TRUE));
      verbose && cat(verbose, "Fields: ", paste(fields, collapse=", "));
      verbose && exit(verbose);
    }
    res <- res[fields];
  }

  verbose && exit(verbose);

  res;
})


############################################################################
# HISTORY:
# 2008-08-12
# o Adding more verbose output.
# 2008-04-13
# o Now readDataFrame() of AffymetrixCdfFile adds "virtual" fields, e.g.
#   the field 'isPm' is inferred and generated from 'pbase' and 'tbase'.
# 2008-04-03
# o Created.
############################################################################
