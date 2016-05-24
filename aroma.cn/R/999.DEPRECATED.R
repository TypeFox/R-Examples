# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# DEFUNCT
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Defunct aroma.cn 1.6.1 (2015-09-18)
setMethodS3("callPeaks", "data.frame", function(fit, ...) {
  ## Defunct 2015-09-18
  .Defunct(msg="callPeaks() for data.frame:s is defunct; use ditto for PeaksAndValleys objects. It will eventually be removed from the package.");

  # Argument 'fit';
  stopifnot(all(is.element(c("type", "x", "density"), colnames(fit))));
  class(fit) <- c("PeaksAndValleys", class(fit));
  callPeaks(fit, ...);
}, private=TRUE, deprecated=TRUE)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# DEPRECATED
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
