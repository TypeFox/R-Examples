# Adds '...' to some base functions. These will later be
# turned into default functions by setMethodS3().

# Needed as of R v 3.1.0 devel (2013-09-21).
nrow <- appendVarArgs(nrow)
ncol <- appendVarArgs(ncol)
flush <- appendVarArgs(flush)
colnames <- appendVarArgs(colnames)
rownames <- appendVarArgs(rownames)

# USED TO DO: rowSums <- appendVarArgs(rowSums)
rowSums <- function(...) UseMethod("rowSums");
setMethodS3("rowSums", "default", function(...) {
  base::rowSums(...);
})

# USED TO DO: rowMeans <- appendVarArgs(rowMeans)
rowMeans <- function(...) UseMethod("rowMeans");
setMethodS3("rowMeans", "default", function(...) {
  base::rowMeans(...);
})

############################################################################
# HISTORY:
# 2006-01-23
# o Created to please R CMD check.
############################################################################
