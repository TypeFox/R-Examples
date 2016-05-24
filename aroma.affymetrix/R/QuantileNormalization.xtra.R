setMethodS3("getOutputIdentifier", "QuantileNormalization", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Calculating the output identifier");

  verbose && enter(verbose, "Retrieving the identifier for input data set");
  ds <- getInputDataSet(this);
  inputId <- getIdentifier(ds);
  verbose && exit(verbose);

  verbose && enter(verbose, "Calculating the identifier for parameters");
  paramId <- this$.paramId;
  params <- getParameters(this);
  params$.targetDistribution <- attr(params$.targetDistribution, "identifier");
  paramId <- getChecksum(list(params));
  this$.paramId <- paramId;
  verbose && exit(verbose);

  verbose && enter(verbose, "Calculating the joint identifier");
  id <- getChecksum(list(inputId, paramId));
  verbose && exit(verbose);

  verbose && exit(verbose);

  id;
}, private=TRUE)



############################################################################
# HISTORY:
# 2006-12-08
# o Should getOutputIdentifier() be made deprecated?
############################################################################
