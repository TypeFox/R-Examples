setConstructorS3("AromaPackage", function(pkgName="aroma.core", ...) {
  extend(Package(pkgName), "AromaPackage");
})


setMethodS3("fixSearchPathInternal", "AromaPackage", function(this, aheadPkgs, behindPkgs, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validating arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # Argument 'aheadPkgs':
  aheadPkgs <- Arguments$getCharacters(aheadPkgs);

  # Argument 'behindPkgs':
  behindPkgs <- Arguments$getCharacters(behindPkgs);


  verbose && enter(verbose, "Reshuffling search path to make it compatible with ", getName(this));

  verbose && cat(verbose, "Packages that should be ahead:");
  verbose && print(verbose, aheadPkgs);
  verbose && cat(verbose, "Packages that should be behind:");
  verbose && print(verbose, behindPkgs);

  # Nothing to do?
  if (length(aheadPkgs) == 0 || length(behindPkgs) == 0) {
    verbose && cat(verbose, "Nothing to do: No packages.");
    verbose && exit(verbose);
    return();
  }

  # Identify the "ahead" package that is last on the search path
  idxs <- match(sprintf("package:%s", aheadPkgs), search());
  lastPkg <- aheadPkgs[which.max(idxs)];
  toPath <- sprintf("package:%s", lastPkg);
  verbose && cat(verbose, "The \"ahead\" package that is last on the search path: ", lastPkg);

  verbose && cat(verbose, "Search path before reshuffling:");
  verbose && print(verbose, search());

  verbose && enter(verbose, "Moving \"behind\" packages");
  # Move those package, if they are loaded.
  pkgsMoved <- c();
  for (pkg in behindPkgs) {
    path <- sprintf("package:%s", pkg);
    if (is.element(path, search())) {
      # Need to move?
      from <- match(path, search());
      to <- match(toPath, search());
      if (from < to) {
        verbose && printf(verbose, "Moving package: %s (%s)\n", pkg, path);
        pkgsMoved <- c(pkgsMoved,
                    moveInSearchPath(from=path, to=toPath, where="after"));
      }
    }
  } # for (pkg ...)
  verbose && exit(verbose);

  verbose && cat(verbose, "Search path after reshuffling:");
  verbose && print(verbose, search());

  verbose && exit(verbose);

  invisible(pkgsMoved);
}, protected=TRUE)

setMethodS3("fixSearchPath", "AromaPackage", function(this, ...) {
})




############################################################################
# HISTORY:
# 2009-05-13
# o Added protected fixSearchPathInternal().
# o Created from 999.AromaAffymetrix.R.
############################################################################
