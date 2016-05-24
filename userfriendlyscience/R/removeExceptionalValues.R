removeExceptionalValues <- function(dat, items=NULL, exception=.005,
                                    silent=FALSE, stringsAsFactors=FALSE) {
  if (is.data.frame(dat)) {
    if (is.null(items)) {
      items <- names(dat);
      if (!silent) {
        cat("No items specified: extracting all variable names in dataframe.\n");
      }
    }
    return(data.frame(lapply(dat, function(x) {
      if (is.numeric(x)) {
        return(ifelse(exceptionalScore(x, prob = exception), NA, x));
      } else {
        return(x);
      }
    }), stringsAsFactors=stringsAsFactors));
  } else {
    return(ifelse(exceptionalScore(dat, prob = exception), NA, dat));
  }
}