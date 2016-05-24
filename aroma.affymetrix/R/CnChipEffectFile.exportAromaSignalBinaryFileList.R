setMethodS3("exportAromaSignalBinaryFileList", "CnChipEffectFile", function(this, whats=c("total", "freqB"), ..., force=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'what':
  whats <- match.arg(whats, several.ok=TRUE);

  if ("freqB" %in% whats) {
    if (!force && this$combineAlleles) {
      whats <- setdiff(whats, "freqB");
      warnings("Ignoring value \"freqB\" in argument 'what', because source does not contain allele-specific chip effects. Use argument 'force=TRUE' to override this and generate a freqB file that contains all NAs: ", getPathname(this));
    }
  }

  NextMethod("exportAromaSignalBinaryFileList", whats=whats);
}, protected=TRUE)


############################################################################
# HISTORY:
# 2008-06-25
# o Created.
############################################################################ 
