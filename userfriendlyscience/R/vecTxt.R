vecTxt <- function(vector, delimiter = ", ", useQuote = "",
                   firstDelimiter = NULL, lastDelimiter = " & ",
                   firstElements = 0, lastElements = 1,
                   lastHasPrecedence = TRUE) {
  
  vector <- paste0(useQuote, vector, useQuote);
  
  if (length(vector) == 1) {
    return(vector);
  }
  
  if (firstElements + lastElements > length(vector)) {
    if (lastHasPrecedence) {
      firstElements <- length(vector) - lastElements;
    } else {
      lastElements <- length(vector) - firstElements;
    }
  }
  
  firstTxt <- lastTxt <- "";
  
  if (is.null(firstDelimiter)) {
    firstDelimiter <- delimiter;
  }
  if (is.null(lastDelimiter)) {
    lastDelimiter <- delimiter;
  }
  
  midBit <- vector;
  if (firstElements > 0) {
    firstBit <- head(vector, firstElements);
    midBit <- tail(vector, -firstElements);
    firstTxt <- paste0(paste0(firstBit, collapse=firstDelimiter), firstDelimiter);
  }
  if (lastElements > 0) {
    lastBit <- tail(vector, lastElements);
    midBit <- head(midBit, -lastElements);
    lastTxt <- paste0(lastDelimiter, paste0(lastBit, collapse=lastDelimiter));
  }
  
  midTxt <- paste0(midBit, collapse=delimiter);
  
  return(paste0(firstTxt, midTxt, lastTxt));
  
}
