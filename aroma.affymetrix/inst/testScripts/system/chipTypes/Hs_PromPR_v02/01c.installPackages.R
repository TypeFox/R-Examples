path <- system.file("testScripts/R", package="aroma.affymetrix");
pathname <- file.path(path, "installUtils.R");
source(pathname);

library("R.utils");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Install
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
verbose && enter(verbose, "Installing test-specific packages");

pkgs <- c("BioC:limma", "gsmoothr");
for (pkg in pkgs) {
  verbose && cat(verbose, "Package: ", pkg);
  installPkg(pkg);
}

verbose && exit(verbose);
