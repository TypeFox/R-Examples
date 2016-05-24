setMethodS3("getProbePositionEffectDesignMatrix", "raw", function(seqs, B=NULL, intercept=TRUE, df=5, ..., verbose=FALSE) {
  # Argument 'seqs':
  if (!is.matrix(seqs)) {
    throw("Argument 'seqs' is of type raw, but not a matrix.");
  }
  P <- ncol(seqs);  # Number of positions in sequences

  # Argument 'df':
  df <- Arguments$getInteger(df, range=c(1,20));

  # Argument 'B':
  if (is.null(B)) {
    # Default: a basis matrix for natural-cubic splines
    B <- splines::ns(1:P, df=df);
  } else if (is.matrix(B)) {
    if (nrow(B) != P) {
      throw("The number of rows in the base-vector matrix does not match the number of positions in the sequences: ", nrow(B), " != ", P);
    }
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup the design matrix
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Building design matrix");

  K <- nrow(seqs);  # Number of sequences
  df <- ncol(B);    # Number of basis vectors
  verbose && cat(verbose, "Number of sequences: ", K);
  verbose && cat(verbose, "Degrees of freedom per position: ", df);

  # Exclude NA factors and last factor
  map <- attr(seqs, "map");
  # Sanity check
  if (is.null(map)) {
    throw("Cannot build design matrix: No attribute 'map'.");
  }
  factors <- map[seq(from=2, to=length(map)-1)];
  verbose && cat(verbose, "Number of factors: ", length(factors));

  L <- df*length(factors);
  if (intercept)
    L <- L + 1;

  # Allocate an KxL prediction matrix in the model
  verbose && enter(verbose, "Allocating design matrix");
  verbose && cat(verbose, "Dimension: ", K, "x", L);
  X <- matrix(0, nrow=K, ncol=L);
  verbose && exit(verbose);

  # Intercept?
  if (intercept) {
    X[,1] <- 1;
    gc <- gc();
  }

  for (bb in seq_along(factors)) {
    verbose && enter(verbose, sprintf("Factor #%d ('%s') of %d",
                               bb, names(factors)[bb], length(factors)));

    # For every position in the sequences
    for (pp in 1:P) {
      verbose && enter(verbose, sprintf("Position #%d of %d", pp, P));

      # Identify sequences with factor 'bb' in position 'pp'
      idxs <- which(seqs[,pp] == factors[bb]);
#      verbose && cat(verbose, "Matching sequences:");
#      verbose && str(verbose, idxs);

      # For every dimension in the base vector
      for (jj in 1:df) {
        cc <- 1 + df*(bb-1) + jj;
        X[idxs,cc] <- X[idxs,cc] + B[pp,jj];
      }
      # Not needed anymore
      idxs <- NULL;

      verbose && exit(verbose);
    } # for (pp ...)

    gc <- gc();
    verbose && exit(verbose);
  } # for (bb ...)

  # Not needed anymore
  seqs <- NULL;
  gc <- gc();

  verbose && cat(verbose, "Design matrix:");
  verbose && str(verbose, X);

  verbose && exit(verbose);

  res <- list(X=X, map=map, factors=factors, B=B);
  class(res) <- "ProbePositionEffectDesignMatrix";

  res;
}) # getProbePositionEffectDesignMatrix()




setMethodS3("getProbePositionEffectDesignMatrix", "character", function(seqs, ..., verbose=FALSE) {
  # Argument 'seqs':
  P <- nchar(seqs);
  P <- unique(P);
  if (length(P) != 1) {
    throw("Argument 'seqs' contains sequences of different lengths: ",
                                   paste(head(sort(P)), collapse=", "));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Coercing sequences into a raw sequence matrix
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Coercing sequences into a raw sequence matrix");
  K <- length(seqs);
  P <- nchar(seqs[1]);
  verbose && cat(verbose, "Number of sequences: ", K);
  verbose && cat(verbose, "object.size(seqs):");
  verbose && print(verbose, object.size(seqs));
  seqs <- paste(seqs, collapse="");
  seqs <- charToRaw(seqs);
  map <- as.raw(0:4);
  names(map) <- c(NA, "A", "C", "G", "T");
  values <- map[2:5];
  from <- charToRaw(paste(names(values), collapse=""));
  for (kk in seq_along(values)) {
    idxs <- which(seqs == from[kk]);
    seqs[idxs] <- values[kk];
  }
  seqs[seqs > map[length(map)]] <- as.raw(0);
  seqs <- matrix(seqs, nrow=K, ncol=P, byrow=TRUE);
  gc <- gc();
  verbose && print(verbose, gc);
  verbose && cat(verbose, "object.size(seqs):");
  verbose && print(verbose, object.size(seqs));
  verbose && exit(verbose);

  getProbePositionEffectDesignMatrix(seqs, ..., verbose=verbose);
})



############################################################################
# HISTORY:
# 2008-07-11
# o Updated to work with the new AromaCellSequenceFile encodings.
# 2008-07-07
# o Added predict().
# o Added getProbePositionEffects().
# o Added na.rm=TRUE to fitProbePositionEffects().
# 2008-07-05
# o Created.
############################################################################
