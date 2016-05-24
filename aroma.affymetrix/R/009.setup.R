.setupAromaAffymetrix <- function(pkg, ...) {
  # To please R CMD check
  ns <- getNamespace("aroma.core");
  .requireBiocPackage <- get(".requireBiocPackage", mode="function", envir=ns);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Patches
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # None at the moment.

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Bioconductor packages aroma.light and affxparser
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # require("aroma.light") - install if missing
  .requireBiocPackage("aroma.light", neededBy=getName(pkg));

  # require("affxparser") - install if missing
  .requireBiocPackage("affxparser", neededBy=getName(pkg));

  # Make sure 'affxparser' is after 'aroma.affymetrix' on the search path
  from <- "package:affxparser";
  to <- "package:aroma.affymetrix";
  fromIdx <- match(from, search());
  toIdx <- match(to, search());
  if (all(is.finite(c(fromIdx, toIdx))) && fromIdx < toIdx) {
    moveInSearchPath(from=from, to=to, where="after");
  }

  # Add custom findCdf() function to affxparser.  This is need to be
  # able to locate CDFs in annotationData/chipTypes/<chipType>/.
  setCustomFindCdf(function(...) {
    AffymetrixCdfFile$findByChipType(..., .useAffxparser=FALSE);
  });


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Package settings (settings might change)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # This code will update the settings according to the default ones.
  updateSettings(pkg);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fix the search path every time a package is loaded
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  setHook("base::library:onLoad", function(...) {
    # Fix the search path
    pkgs <- fixSearchPath(aroma.affymetrix);
    if (length(pkgs) > 0) {
      warning("Packages reordered in search path: ",
                                            paste(pkgs, collapse=", "));
    }
  }, action="append");
} # .setupAromaAffymetrix()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# affxparser related
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Instead of asking users to write affxparser::writeCdf() when
# aroma.affymetrix is loaded...
writeCdf.default <- function(...) {
  ns <- loadNamespace("affxparser");
  `affxparser::writeCdf` <- get("writeCdf", envir=ns, mode="function");
  `affxparser::writeCdf`(...)
}



############################################################################
# HISTORY:
# 2013-08-03
# o Added writeCdf.default() which loads namespace 'affxparser' and calls
#   writeCdf() of affxparser, i.e. effectively affxparser::writeCdf().
# 2009-02-22
# o Removed code checking for package updates. Was never called.
# 2009-02-21
# o Added setting memory$ram=1.
# 2008-02-14
# o Renamed existing threshold hold to 'timestampsThreshold',
#   'medianPolishThreshold', and 'skipThreshold'.
# 2008-02-12
# o Added default values for settings 'models$RmaPlm$...'.
# 2008-01-30
# o Added default values for settings 'rules$allowAsciiCdfs' and
#   'output$maxNbrOfArraysForTimestamps'.
# 2007-12-13
# o Added code for automatic updates on startup.  In active by default.
# o Added settings for 'checkForPatches' and 'checkInterval'.
# o Now the settings are set according to a tempate, if missing.
# 2007-08-30
# o Added "patch" to make sure that there is rowMedians() supporting
#   missing values.
# 2007-07-04
# o Removed the patch for digest(); digest v0.3.0 solved the problem.
# o Added a patch of functions in 'base', e.g. matrix().
# 2007-04-04
# o Moved the patch of digest() here.
# 2007-03-07
# o Added test for consistency of digest().
# 2007-03-06
# o Added onLoad hook function for library() and require() to call
#   fixSearchPath() of the package, which reorders the search path so that
#   problematic packages are after this package in the search path.
# 2007-02-22
# o Added default settings stubs.
# 2007-02-12
# o Created.
############################################################################
