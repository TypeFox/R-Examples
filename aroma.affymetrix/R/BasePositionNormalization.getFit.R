setMethodS3("getFit", "BasePositionNormalization", function(this, array, ...) {
  # Assert that the normalization has been done
  if (!isDone(this)) {
    throw("Cannot get fit: Data is not processed: ", getFullName(this));
  }

  # Argument 'array':
  inSet <- getInputDataSet(this);
  array <- Arguments$getIndices(array, max=length(inSet));

  outSet <- getOutputDataSet(this);
  outDf <- outSet[[array]];
  path <- getPath(this);
  filename <- sprintf("%s,fit.RData", getFullName(outDf));
  pathname <- Arguments$getReadablePathname(filename, path=path, mustExist=TRUE);
  res <- loadObject(pathname);
  fit <- res$fit;
  # Not needed anymore
  res <- NULL;

  fit;
}, protected=TRUE);



############################################################################
# HISTORY:
# 2008-12-15
# o Add getFit().
############################################################################
