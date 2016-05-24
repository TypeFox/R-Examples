setMethodS3("setAttributesBy", "AromaMicroarrayDataSet", function(this, object, ...) {
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


setMethodS3("setAttributesBySampleAnnotationSet", "AromaMicroarrayDataSet", function(this, sas, ..., verbose=FALSE) {
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


setMethodS3("setAttributesBySampleAnnotationFile", "AromaMicroarrayDataSet", function(this, saf, force=FALSE, ..., verbose=FALSE) {
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
  res <- applyTo(saf, names, FUN=setAttrs, force=force, verbose=verbose);

  invisible(this);
}, protected=TRUE)


############################################################################
# HISTORY:
# o Moved attribute methods from AffymetrixCelSet to AromaMicroarrayDataSet.
# 2007-03-24
# o Now the setAttributeByNnn() methods return itself invisibly.
# 2007-03-14
# o Now setAttributesBySampleAnnotationFile() also set attributes.
# 2007-03-06
# o Added setAttributesBy().
# o Added setAttributesBySampleAnnotationFile().
# o Created.
############################################################################
