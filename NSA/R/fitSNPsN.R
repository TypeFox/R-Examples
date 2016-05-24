fitSNPsN <- function(T, references, threshold =.135, ..., verbose = verbose) {
  # This is an internal function. Because of this, we will assume that
  # all arguments are valid and correct.  No validation will be done.

  thr <- references > threshold;

  # check if there are SNPs with no references, then use all of them
  aux <- rowSums(!thr);
  thr[aux==ncol(thr),] <- TRUE;
  
  ind <- thr == FALSE;
  thr[ind] <- NA;
  
  #if(sum(is.na(thr))==length(thr)){
  #  thr[1:length(thr)] <- TRUE;
  #}

  ref <- rowMedians(T*thr, na.rm=TRUE);

  res <- 2*T/ref;
  res;
} # fitSNPSN()

###########################################################################
# HISTORY:
# 2010-06-24 [MO]
# o Created from fitCalMaTe.R.
###########################################################################
