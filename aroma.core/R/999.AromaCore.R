setConstructorS3("AromaCore", function(pkgName="aroma.core", ...) {
  extend(AromaPackage(pkgName), "AromaCore");
})




setMethodS3("fixSearchPath", "AromaCore", function(this, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # RULES
  # 2009-05-13:
  # o grid must be after aroma.core. [getNames()]
  # 2008-08-27:
  # o affy must be after aroma.light, otherwise the former overrides
  #   the generic plotDensity() function of the latter.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Figure out which of our packages (aroma.core, aroma.light etc.) is
  # last on the search path.
  aheadPkgs <- c("aroma.core", "aroma.light");

  # Problematic package that must be after this package on the search path
  behindPkgs <- c("grid", "affy");

  res <- fixSearchPathInternal(this, aheadPkgs=aheadPkgs,
                                           behindPkgs=behindPkgs, ...);

  # Return the package actually moved
  invisible(res);
})



############################################################################
# HISTORY:
# 2009-05-13
# o Now extending the AromaPackage class.
# o Created from 999.AromaAffymetrix.R.
############################################################################
