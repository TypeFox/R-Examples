path <- system.file("testScripts/R", package="aroma.affymetrix");
pathname <- file.path(path, "installUtils.R");
source(pathname);

library("R.utils");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Install
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
verbose && enter(verbose, "Installing test-specific packages");

pkgs <- c("BioC:gcrma")
pkgs <- c(pkgs, "BioC:affy", "BioC:limma");
pkgs <- c(pkgs, "BioC:oligo");
pkgs <- c(pkgs, "BioC:hgu133plus2cdf");
pkgs <- c(pkgs, "BioC:hgu133plus2.db");
pkgs <- c(pkgs, "BioC:pd.hg.u133.plus.2");

for (pkg in pkgs) {
  verbose && cat(verbose, "Package: ", pkg);
  installPkg(pkg);
}

verbose && exit(verbose);
