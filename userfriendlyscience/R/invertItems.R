invertItems <- function(dat, items = NULL, ...) {
  if (is.null(items)) {
    items <- names(dat);
  } else if ((!is.character(items)) && (!is.numeric(items))) {
    stop("Argument 'items' is not a character vector or numeric vector ",
         "(but instead of type ", typeof(dat), ").");
  }
  usedDat <- dat[, items];
  
  ### Previous inversions
  prevInv <- lapply(dat[, items], attr, 'inverted');
  ### Replace NULL with FALSE
  prevInv <- lapply(prevInv, function(x) ifelse(is.null(x), FALSE, x));
  ### Warn if one or more items were already inverted
  if (sum(unlist(prevInv)) > 0) {
    alreadyInverted <- names(prevInv)[unlist(prevInv)];
    warning("Variables (columns) ", vecTxt(alreadyInverted, useQuote='"'),
            " have already been inverted! ",
            "Set ignorePreviousInversion to TRUE to override this check ",
            "and invert the vector anyway.");
    usedDat <- usedDat[, !(names(usedDat) %in% alreadyInverted)];
  }
  
  ### All convert factors to numeric vectors
  usedDat <- massConvertToNumeric(usedDat);
  
  ### Check whether any non-numeric vectors remain
  invalidVectors <- lapply(usedDat, is.numeric);
  if (FALSE %in% unlist(invalidVectors)) {
    invalidVectors <- names(invalidVectors)[!unlist(invalidVectors)];
    warning("Variables (columsn) ", vecTxt(invalidVectors, useQuote='"'),
            " have a type other than numeric or factor! Ignoring these.");
    usedDat <- usedDat[, !(names(usedDat) %in% invalidVectors)];
  }
  
  items <- names(usedDat);
  
  dat[, items] <- data.frame(lapply(usedDat, invertItem), ...);
  return(dat);
}