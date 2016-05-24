### Convert a vector to numeric values and trying to be smart about it.
convertToNumeric <- function (vector, byFactorLabel = FALSE) {
  ### Check whether the vector is datetime
  if (sum(sapply(class(vector), grepl, pattern='POSIX')) > 0) {
    return(vector);
  }
  if (!(is.factor(vector) | is.numeric(vector) |
        is.character(vector) | is.logical(vector))) {
    stop("Argument 'vector' must be a vector! Current class = '",
         class(vector), "'. To mass convert e.g. a dataframe, ",
         "use massConvertToNumber.");
  }
  if(is.factor(vector) && byFactorLabel) {
    ### Decimal symbol might be a comma instead of a period: convert
    ### factor to character vector and replace commas with periods
    vector <- as.numeric(gsub(as.character(vector), pattern=",", replacement="."));
    return();
  }
  else if (is.character(vector)) {
    return(suppressWarnings(as.numeric(gsub(as.character(vector), pattern=",", replacement="."))));
  }
  else {
    ### Thus, for numeric vectors; factors to be converted by index of the levels
    ### instead of by their labels; and logical vectors.
    return(as.numeric(vector));
  }
}
