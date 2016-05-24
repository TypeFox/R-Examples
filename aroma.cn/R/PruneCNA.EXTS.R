setMethodS3("nbrOfGenerations", "PruneCNA", function(this, ...) {
  length(this);
})

setMethodS3("extractGenerations", "PruneCNA", function(this, generations, ...) {
  # Argument 'generations':
  n <- nbrOfGenerations(this);
  if (any(generations < 0)) {
    generations <- setdiff(seq_len(n), -generations);
  }
  generations <- Arguments$getIndices(generations, max=n);

  res <- unclass(this);
  res <- res[generations];
  class(res) <- class(this);
  res;
})

setMethodS3("[", "PruneCNA", function(x, i) {
  # To please R CMD check
  this <- x;

  extractGenerations(this, generations=i);
})


############################################################################
# HISTORY:
# 2011-02-06
# o Added nbrOfGenerations(), extractGenerations() and "[" for PruneCNA.
# o BUG FIX: The test for identifying atomic islands from a set of
#   atomic regions assumed all atomic regions where of length H=1.
############################################################################
