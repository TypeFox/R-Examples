installPkg <- function(pkg, ...) {
  # Parse 'pkg'
  pattern <- "^([^:]*):([^:]*)$";
  if (regexpr(pattern, pkg) != -1) {
    repos <- gsub(pattern, "\\1", pkg);
    pkg <- gsub(pattern, "\\2", pkg);
  } else {
    repos <- "CRAN";
  }

  # Nothing to do?
  if (isPackageInstalled(pkg)) return();

  # Assert known repository
  repos <- match.arg(repos, choices=c("CRAN", "BioC", "R-Forge"));

  tryCatch({
    if (repos == "CRAN") {
      install.packages(pkg);
    } else if (repos == "BioC") {
      source("http://www.bioconductor.org/biocLite.R");
      biocLite(pkg, suppressUpdates=TRUE, ask=FALSE);
    } else if (repos == "R-Forge") {
      install.packages(pkg, repos="http://R-Forge.R-project.org");
    }
  }, error = function(ex) {
    print(ex);
  })
} # installPkg()


############################################################################
# HISTORY:
# 2012-09-01
# o Added installPkg(), which installs *some* known packages by their
#   names without having to specify the repository.
# o Created.
############################################################################
