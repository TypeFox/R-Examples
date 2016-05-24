# Truncate ASCNs to avoid non-positives.
# Truncation is done such that TCN is preserved regardlessly.
setMethodS3("truncateThetaAB", "array", function(data, ...) {
  # This is an internal function. Because of this, we will assume that
  # all arguments are valid and correct.  No validation will be done.
  dim <- dim(data);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Estimate corrections
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (CA*,CB*) = (CA,CB) + (-dA,+dA) for all SNPs where CA < 0.
  x <- data[,1,];
  idxsA <- which(x < 0);
  dA <- x[idxsA];
  rm(x);

  # (CA*,CB*) = (CA,CB) + (+dB,-dB) for all SNPs where CB < 0.
  x <- data[,2,];
  idxsB <- which(x < 0);
  dB <- x[idxsB];
  rm(x);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Apply coorections
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  delta <- array(0, dim=dim[-2]);
  delta[idxsA] <- dA;
  rm(idxsA,dA);
  data[,1,] <- data[,1,] - delta;
  data[,2,] <- data[,2,] + delta;
  rm(delta);

  delta <- array(0, dim=dim[-2]);
  delta[idxsB] <- dB;
  rm(idxsB,dB);
  data[,1,] <- data[,1,] + delta;
  data[,2,] <- data[,2,] - delta;
  rm(delta);

  data;
}) # truncateThetaAB()



setMethodS3("truncateThetaAB", "matrix", function(data, ...) {
  # This is an internal function. Because of this, we will assume that
  # all arguments are valid and correct.  No validation will be done.

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Estimate corrections
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (CA*,CB*) = (CA,CB) + (-dA,+dA) for all SNPs where CA < 0.
  x <- data[1,];
  idxsA <- which(x < 0);
  dA <- x[idxsA];
  rm(x);

  # (CA*,CB*) = (CA,CB) + (+dB,-dB) for all SNPs where CB < 0.
  x <- data[2,];
  idxsB <- which(x < 0);
  dB <- x[idxsB];
  rm(x);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Apply coorections
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  delta <- array(0, dim=ncol(data));
  delta[idxsA] <- dA;
  rm(idxsA,dA);
  data[1,] <- data[1,] - delta;
  data[2,] <- data[2,] + delta;
  rm(delta);

  delta <- array(0, dim=ncol(data));
  delta[idxsB] <- dB;
  rm(idxsB,dB);
  data[1,] <- data[1,] + delta;
  data[2,] <- data[2,] - delta;
  rm(delta);

  data;
}) # truncateThetaAB()


###########################################################################
# HISTORY:
# 2010-06-22 [MO]
# o Added a truncateThetaAB() for matrix.
# 2010-06-19 [HB]
# o Created.
###########################################################################
