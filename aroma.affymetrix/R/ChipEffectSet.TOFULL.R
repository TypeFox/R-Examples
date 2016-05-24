setMethodS3("getAsFullCelSet", "ChipEffectSet", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Getting chip effect set expanded to the full CDF");
  cells <- NULL;
  files <- list();
  for (kk in seq_along(this)) {
    cef <- this[[kk]];
    verbose && enter(verbose, sprintf("Array #%d ('%s') %d",
                                          kk, getName(cef), length(this)));

    cf <- getAsFullCelFile(cef, ..., cells=cells, verbose=less(verbose, 5));

    cells <- attr(cf, "cells");
    attr(cf, "cells") <- NULL;

    files[[kk]] <- cf;
    verbose && exit(verbose);
  }
  verbose && exit(verbose);

  res <- AffymetrixCelSet(files);

  res;
}, protected=TRUE);



############################################################################
# HISTORY:
# 2008-03-31
# o Now the cell map is reused across arrays.
# 2008-03-20
# o Created.
############################################################################
