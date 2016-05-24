setCustomFindCdf <- function(fcn, ...) {
  # Argument 'fcn':
  if (!is.function(fcn)) {
    stop("Argument 'fcn' is not a function: ", mode(fcn));
  }

  settings <- getOption("affxparser.settings");
  if (is.null(settings))
    settings <- list();

  oldFcn <- settings$methods$findCdf;
  settings$methods$findCdf <- fcn;
  options("affxparser.settings"=settings);

  invisible(oldFcn);
} # setCustomFindCdf()


############################################################################
# HISTORY:
# 2007-02-12
# o Created. To be added to affxparser later.
############################################################################  
