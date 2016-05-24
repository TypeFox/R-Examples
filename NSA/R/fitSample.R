fitSample <- function(T, references, threshold =.135, ..., verbose = verbose) {
  # This is an internal function. Because of this, we will assume that
  # all arguments are valid and correct.  No validation will be done.

  thr <- references > threshold;
  thr[thr==FALSE] <- NA;
  
  if(sum(is.na(thr))==length(thr)){
    thr[1:length(thr)] <- TRUE;
  }
  
  ref <- median(T*thr, na.rm=TRUE);
  res <- 2*T/ref;
  
  res;
} # fitSample()

###########################################################################
# HISTORY:
# 2010-06-29 [MO]
# o Created from fitSNPSN.R.
###########################################################################
