setMethodS3("getChromosomes", "CopyNumberRegions", function(this, ...) {
  sort(unique(this$chromosome));
})


setMethodS3("nbrOfChromosomes", "CopyNumberRegions", function(this, ...) {
  length(getChromosomes(this));
})


setMethodS3("unique", "CopyNumberRegions", function(x, ...) {
  # To please R CMD check
  this <- x;

  # Keep only unique regions
  data <- as.data.frame(this);
  data <- unique(data);  
  fields <- colnames(data);
  res <- clone(this);
  for (field in fields) {
    res[[field]] <- data[[field]];
  }

  res;
})

setMethodS3("append", "CopyNumberRegions", function(this, other, ...) {
  # Argument 'other':
  other <- Arguments$getInstanceOf(other, class(this)[1]);

  data <- as.data.frame(this);
  fields <- colnames(data);
  for (field in fields) {
    values <- c(this[[field]], other[[field]]);
    this[[field]] <- values;
  }

  invisible(this);
})

############################################################################
# HISTORY:
# 2010-08-05
# o Created.
############################################################################
