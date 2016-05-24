setMethodS3("getChromosomes", "Arguments", function(static, chromosomes, range=c(1,24), ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'chromosome':
  if (is.character(chromosomes)) {
    # Interpret X and Y
    map <- c("X"=23, "Y"=24);
    idxs <- (chromosomes %in% names(map));
    chromosomes[idxs] <- map[chromosomes[idxs]];
  }

  chromosomes <- Arguments$getIndices(chromosomes, range=range, ...);

  chromosomes;
}, static=TRUE, protected=TRUE)


setMethodS3("getChromosome", "Arguments", function(static, chromosome, ..., length=1) {
  getChromosomes(static, chromosomes=chromosome, length=length, ...);
}, static=TRUE, protected=TRUE)


############################################################################
# HISTORY:
# 2011-02-18
# o Moved static getTags() to Arguments to R.filesets v0.9.3.
# 2010-01-25
# o Added static getTags() to Arguments.
# 2007-03-15
# o BUG FIX: Forgot to pass '...' down in getChromosome().
# 2007-03-06
# o Added getChromosome() and getChromosomes().
# o Created.
############################################################################
