readCdfGroupStrands <- function(..., what=c("probe", "target")) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'what':
  what <- match.arg(what);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read the strand information from the CDF file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- .readCdf(..., readUnitDirection=TRUE, readGroupDirection=TRUE,
         readXY=FALSE, readBases=FALSE, readIndexpos=FALSE, readIsPm=FALSE,
                                      readAtoms=FALSE, readUnitType=FALSE);

  getTarget <- (what == "target");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract the strandinformation
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  swapIdxs <- c();
  for (uu in seq_along(cdf)) {
    unit <- .subset2(cdf, uu);

    unitDirection <- .subset2(unit, "unitdirection");
    strands <- unlist(.subset2(unit, "groups"), use.names=FALSE);
    if (identical(unitDirection, "antisense")) {
      strands <- c(sense="antisense", antisense="sense")[strands];
      swapIdxs <- c(swapIdxs, uu);
    }

    # Get strands for target?
    if (getTarget)
      strands <- c(sense="antisense", antisense="sense")[strands];
    strands <- unname(strands);

    cdf[[uu]] <- strands;
  }

  attr(cdf, "swapIdxs") <- swapIdxs;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Return a list of character vector
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf;
}

############################################################################
# HISTORY:
# 2006-12-12
# o Created. Should be added to affxparser.
############################################################################
