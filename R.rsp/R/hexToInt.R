
hexToInt <- function(hex, ...) {
  hex <- as.character(hex);
  hex <- tolower(hex);
  hex <- strsplit(hex, split="", fixed=TRUE)[[1]];
  hex <- match(hex, c(0:9, letters[1:6]))-1;
  sum(16^seq(from=length(hex)-1,to=0,by=-1) * hex);
}


###############################################################################
# HISTORY:
# 2005-09-24
# o Created.
###############################################################################

