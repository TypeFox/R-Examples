### Function to format p values nicely
formatPvalue <- function (values, digits = 3, spaces=TRUE, includeP = TRUE) {
  missingValues <- is.na(values);
  values <- ifelse(values < 0, 0, ifelse(values > 1, 1, values));
  pchar <- ifelse(includeP, "p = ", "");
  eps <- 10 ^ -digits;
  res <- paste0(pchar, noZero(format.pval(round(values, digits),
                                          eps=eps, digits=digits,
                                          scientific=digits+1)));
  if (spaces) {
    res <- gsub("= <", "< ", res);
  } else {
    res <- gsub("= <", "<", res);
    res <- gsub(" ", "", res);
  }
  res <- ifelse(missingValues, NA, res);
  return(res);
}
