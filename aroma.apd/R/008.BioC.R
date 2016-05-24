# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Bioconductor related
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
.requireBiocPackage <- function(package, neededBy="aroma.core", ...) {
  # Trick 'R CMD check' to not generate NOTEs.
  requireX <- base::require;
  catX <- base::cat;

  res <- suppressWarnings({
    requireX(package, character.only=TRUE);
  });

  # Not installed?
  if (!res) {
    if (interactive()) {
      # Trick 'R CMD check' to not generate NOTEs.
      catX("Package '", package, "' is not available or could not be loaded. Will now try to install it from Bioconductor (requires working internet connection):\n");

      # To please R CMD check
      biocLite <- NULL; rm(list="biocLite");
      source("http://www.bioconductor.org/biocLite.R");
      biocLite(package);
      # Assert that the package can be successfully loaded
      res <- requireX(package, character.only=TRUE);
      if (!res) {
        throw("Package 'affxparser' could not be loaded. Please install it from Bioconductor, cf. http://www.bioconductor.org/");
      }
    } else {
      warning("Package '", package, "' could not be loaded. Without it ", neededBy, " will not work. Please install it from Bioconductor, cf. http://www.bioconductor.org/");
    }
  }
} # .requireBiocPackage()


############################################################################
# HISTORY:
# 2013-01-05
# o Copied .requireBiocPackage() from aroma.core.
############################################################################
