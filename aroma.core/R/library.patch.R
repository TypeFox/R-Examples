library <- function(...) {
  res <- withVisible(base::library(...));
  callHooks("base::library:onLoad");
  value <- res$value;
  if (res$visible) {
    return(value);
  } else {
    return(invisible(value));
  }
} # library()

require <- function(...) {
  res <- withVisible(base::require(...));
  callHooks("base::library:onLoad");
  value <- res$value;
  if (res$visible) {
    return(value);
  } else {
    return(invisible(value));
  }
} # library()


############################################################################
# HISTORY:
# 2011-12-22
# o BUG FIX: The overridden library() would always return an invisible()
#   object, even if base::library() wouldn't.
# 2007-03-06
# o Added onLoad hooks to library() and require().
# o Created.
############################################################################

