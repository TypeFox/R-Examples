# Truncate Fracs to avoid negatives and over 1.
# Truncation is done such that TCN is preserved regardlessly.
setMethodS3("truncateFracB", "matrix", function(data, ...) {
  # This is an internal function. Because of this, we will assume that
  # all arguments are valid and correct.  No validation will be done.

  # Truncate fracB values to 0 and 1.
  C <- colSums(data);
  fracB <-  data[2,] / C;
  eps <- 0;
  
  fracB[(fracB < eps)] <- eps;
  fracB[(fracB > 1)] <- 1;    

  data[1,] <- C*(1-fracB);
  data[2,] <- C*(fracB);

  data;
}) # truncateFracB()


# Truncate ASCNs to avoid non-positives.
# Truncation is done such that TCN is preserved regardlessly.
setMethodS3("truncateFracB", "array", function(data, ...) {
  # This is an internal function. Because of this, we will assume that
  # all arguments are valid and correct.  No validation will be done.

  # Truncate fracB values to 0 and 1.
  fracB <-  data[,2,];
  eps <- 0;
  
  fracB[(fracB < eps)] <- eps;
  fracB[(fracB > 1)] <- 1;    

  dataA <- data[,1,]*(1-fracB);
  dataB <- data[,1,]*(fracB);
  
  data[,1,] <- dataA+dataB;
  data[,2,] <- fracB;
  
  data;
}) # truncateFracB()
###########################################################################
# HISTORY:
# 2010-06-22 [MO]
# o Created.
###########################################################################
