setConstructorS3("GenericDataFileList", function(dfList=list(), ...) {
  # Argument 'dfList':
  if (is.list(dfList)) {
    className <- "GenericDataFile";
    for (kk in seq_along(dfList)) {
      df <- dfList[[kk]];
      if (!inherits(df, className)) {
        throw(sprintf("Element #%d of argument 'dfList' is not of class %s: ", kk, class(df)[1]));
      }
    }
  } else {
    throw("Argument 'dfList' is not a list: ", dfList);
  }

  extend(dfList, "GenericDataFileList");
})



###########################################################################
# HISTORY:
# 2009-05-12
# o Created.
###########################################################################
