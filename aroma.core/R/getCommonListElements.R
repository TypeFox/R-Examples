getCommonListElements <- function(lst, ignoreEmpty=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 1. Scan list for common elements
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Vector to store common elements
  common <- c();
  n <- length(lst);
  first <- TRUE;
  for (kk in seq_len(n)) {
    value <- lst[[kk]];

    # Ignoring empty list elements or not?
    if (is.null(value)) {
      if (ignoreEmpty)
        next;

      # We know for sure that the rest will be empty too
      common <- c();
      break;
    }

    if (first) {
      common <- value;
      first <- FALSE;
    } else {
      common <- intersect(common, value);

      # Done?
      if (length(common) == 0)
        break;
    }
  } # for (kk ...)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 2. Keep only common elements
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  lst <- lapply(lst, FUN=function(ss) {
    keep <- (ss %in% common);
    ss[keep];
  });

  lst;
} # getCommonListElements()


############################################################################
# HISTORY:
# 2007-09-26
# o BUG FIX: getCommonListElements() would exclude duplicated tags, e.g. -X.
# 2007-02-16
# o Created.
############################################################################
