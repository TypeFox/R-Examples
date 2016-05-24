massConvertToNumeric <- function (dat, byFactorLabel = FALSE,
                                  ignoreCharacter = TRUE,
                                  stringsAsFactors = FALSE) {
  storedAttributes <- attributes(dat);
  dat <- data.frame(lapply(dat, function(x) {
    if (is.character(x) && ignoreCharacter) {
      return(x);
    }
    else {
      return(convertToNumeric(x, byFactorLabel = byFactorLabel));
    }
  }), stringsAsFactors=stringsAsFactors);
  attributes(dat) <- storedAttributes;
  return(dat);
}