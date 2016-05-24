# Removes all segment calls and corresponding parameter estimates.
setMethodS3("clearCalls", "AbstractCBS", function(fit, ...) {
  segs <- fit$output;
  params <- fit$params;

  # Drop all calls
  excl <- grep("Call$", colnames(segs));
  if (length(excl) > 0L) {
    segs <- segs[,-excl];
  }

  # Drop all call parameters (AD HOC!)
  for (ff in c("deltaROH", "deltaAB", "deltaLOH")) {
    params[[ff]] <- NULL;
  }

  fit$output <- segs;
  fit$params <- params;

  invisible(fit);
}, protected=TRUE)



##############################################################################
# HISTORY
# 2013-10-24
# o Added clearCalls() for AbstractCBS.
# o Created.
##############################################################################
