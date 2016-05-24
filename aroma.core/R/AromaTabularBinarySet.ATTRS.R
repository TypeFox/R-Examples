setMethodS3("setAttributesBy", "AromaTabularBinarySet", function(this, object, ...) {
  methodName <- sprintf("setAttributesBy%s", class(object)[1]);
  if (!exists(methodName, mode="function")) {
    throw("No set function found: ", methodName);
  }

  fcn <- get(methodName, mode="function");
  tryCatch({
    fcn(this, object, ...);
  }, error = function(ex) {
    print(ex);
    throw("Failed to apply attributes by object of class: ", class(object)[1]);
  })

  invisible(this);
}, protected=TRUE)


setMethodS3("setAttributesByTags", "AromaTabularBinarySet", function(this, tags=getTags(this), ...) {
  # Split tags
  tags <- Arguments$getTags(tags, collapse=NULL);

  newAttrs <- NextMethod("setAttributesByTags", tags=tags);

  # Parse XY, XX, XXX etc tags
  values <- grep("^X*Y*$", tags, value=TRUE);
  if (length(values) > 0) {
    newAttrs <- c(newAttrs, setAttributeXY(this, values));
  }

  # Parse tri<chromosome> tags
  values <- grep("^tri([1-9]|[0-9][0-9]|X|Y)$", tags, value=TRUE);
  if (length(values) > 0) {
    values <- gsub("^tri", "", values);
    chromosomes <- Arguments$getChromosomes(values);
    keys <- sprintf("n%02d", chromosomes);
    newAttrs <- c(newAttrs, lapply(keys, FUN=function(key) {
      setAttribute(this, key, 3);
    }));
  }

  # Return nothing
  invisible(newAttrs);
}, protected=TRUE)


setMethodS3("setAttributesBySampleAnnotationSet", "AromaTabularBinarySet", function(this, sas, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  res <- lapply(sas, FUN=function(saf) {
    verbose && enter(verbose, "Applying sample annotations");
    on.exit({verbose && exit(verbose)});

    setAttributesBy(this, saf, ..., verbose=less(verbose));
  });

  invisible(this);
}, protected=TRUE)


setMethodS3("setAttributesBySampleAnnotationFile", "AromaTabularBinarySet", function(this, saf, force=FALSE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  setAttrs <- function(appliesTo, tags=NULL, ..., verbose=FALSE) {
    verbose && enter(verbose, "Applying sample annotations");
    on.exit({verbose && exit(verbose)});

    args <- list(...);
    nargs <- length(args);

    tags <- Arguments$getTags(tags, collapse=NULL);

    if (!is.null(tags)) {
      verbose && cat(verbose, "Tags: ", paste(tags, collapse=", "));
      nargs <- nargs + 1;
    }

    # Nothing to do?
    if (nargs == 0)
      return();

    # Typically the below only applies to one sample
    verbose && cat(verbose, "Applies to ", length(appliesTo), " sample(s).");
    for (kk in seq_along(appliesTo)) {
      idx <- appliesTo[kk];
      verbose && cat(verbose, "Sample: ", names(appliesTo)[kk]);

      # Get the CEL file
      cf <- this[[idx]];

      # Apply the attributes
      setAttributes(cf, ...);

      # Apply the tags
      setAttributesByTags(cf, tags)
    };
  } # setAttrs()

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'saf':
  saf <- Arguments$getInstanceOf(saf, "SampleAnnotationFile");

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  names <- getFullNames(this);
  # AD HOC
  names <- gsub(",total", "", names, fixed=TRUE);
  res <- applyTo(saf, names, FUN=setAttrs, force=force, verbose=verbose);

  invisible(this);
}, protected=TRUE)


############################################################################
# HISTORY:
# 2010-02-17
# o Created.
############################################################################
