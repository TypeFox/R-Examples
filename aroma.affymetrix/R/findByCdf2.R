setMethodS3("findByCdf2", "default", function(chipType, tags=NULL, nbrOfUnits=NULL, validator=NULL, firstOnly=TRUE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'chipType':
  chipType <- Arguments$getCharacter(chipType, length=c(1,1));

  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- unlist(strsplit(tags, split=",", fixed=FALSE), use.names=FALSE);
    tags <- trim(tags);
    tags <- tags[nzchar(tags)];
    if (length(tags) == 0)
      tags <- NULL;
  }

  # Argument 'validator':
  if (!is.null(validator)) {
    if (is.function(validator)) {
    } else {
      throw("Argument 'validator' is not a function: ", mode(validator));
    }
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Generate all possible fullname 'chipTypes' and search for the existance
  # of a CDF with the longest name.
  verbose && enter(verbose, "Searching for CDFs");
  verbose && cat(verbose, "Chip type: ", chipType);
  verbose && cat(verbose, "Tags: ", paste(tags, collapse=", "));

  pathnames <- NULL;

  cdf <- NULL;
  for (kk in rev(c(0,seq_along(tags)))) {
    cdfTags <- tags[seq_len(kk)];
    fullname <- paste(c(chipType, cdfTags), collapse=",");
    verbose && printf(verbose, "Trying '%s'...", fullname);
    tryCatch({
      cdf <- AffymetrixCdfFile$byChipType(chipType, tags=cdfTags, nbrOfUnits=nbrOfUnits);

      # Validate?  If invalid, skip it.
      if (!is.null(validator)) {
        if (!validator(cdf)) {
          cdf <- NULL;
        }
      }
    }, error = function(ex) {})

    # Found a CDF?
    if (!is.null(cdf)) {
      pathnames <- c(pathnames, getPathname(cdf));
      verbose && writeRaw(verbose, "found.\n");
      if (firstOnly)
        break;
    } else {
      verbose && writeRaw(verbose, "no match.\n");
    }
  }
  verbose && exit(verbose);

  pathnames;
}) # findByCdf2()


############################################################################
# HISTORY:
# 2009-02-11
# o Added argument 'nbrOfUnits' to fromCdf2().
# 2008-01-19
# o Created from fromCdf() in AromaUnitTabularBinaryFile.R.
############################################################################
