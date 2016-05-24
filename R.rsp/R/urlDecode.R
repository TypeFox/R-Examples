urlDecode <- function(url, ...) {
  pattern <- "(.*)(%)([[:xdigit:]]{2})(.*)";
  value <- url;
  if (is.na(value))
    return("");
  value <- gsub("+", " ", value, fixed=TRUE);
  while(TRUE) {
    pos <- regexpr(pattern, value);
    if (pos == -1)
      break;
    hex <- sub(pattern, "\\3", value);
    ascii <- intToChar(hexToInt(hex));
    value <- sub(pattern, paste("\\1", ascii, "\\4", sep=""), value);
  }
  value;
}

###############################################################################
# HISTORY:
# 2006-02-22
# o BUG FIX: urlDecode(NA) gave an error. Now it returns "".
# 2005-09-24
# o Created.
###############################################################################
