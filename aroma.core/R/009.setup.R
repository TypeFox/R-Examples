# Dummy to please R CMD check
aromaSettings <- NULL

.loadSettings <- function(pkgname, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read settings file ".<name>Settings" and store it in package
  # variable '<name>Settings'.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # ...but don't load settings if running R CMD check
  name <- "aroma"
  varName <- sprintf("%sSettings", name)
  basename <- paste(".", varName, sep="")
  if (queryRCmdCheck() == "notRunning") {
    settings <- AromaSettings$loadAnywhere(basename, verbose=TRUE)
  } else {
    settings <- NULL
  }
  if (is.null(settings)) {
    settings <- AromaSettings(basename)
  }

  assign(varName, settings, envir=getNamespace(pkgname))
} # .loadSettings()


.setupAromaCore <- function(pkg, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Assert that digest() gives a consistent result across R versions
  # and platforms.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!identical(getOption("aroma.core::assertDigest"), FALSE)) {
    R.cache:::.assertDigest("error")
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fix the search path every time a package is loaded
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  setHook("base::library:onLoad", function(...) {
    # Fix the search path
    pkgs <- fixSearchPath(aroma.core)
    if (length(pkgs) > 0L) {
      warning("Packages reordered in search path: ",
                                            paste(pkgs, collapse=", "))
    }
  }, action="append")
} # .setupAromaCore()



############################################################################
# HISTORY:
# 2012-08-26
# o Added .loadSettings() extracted from .setupAromaCore().  It now
#   assigned the 'aromaSettings' variable when package is loaded.
#   Previously it had to be attach.
# 2012-04-22
# o CLEANUP: Package no longer try to apply package patches, which was
#   only possible when namespaces where not used.
# 2012-01-12
# o CLEANUP: Dropped internal patch of base::serialize(), because it
#   was only applied to R (< 2.12.0) anyway and this package now
#   requires R (>= 2.12.0).
# 2012-01-11
# o ROBUSTNESS: Aroma settings are no longer loaded during R CMD check.
# 2010-10-27
# o CLEANUP: Removed outdated patch for finding smoothScatter(), which
#   is not needed in R v2.9.0 and beyond.
# o CLEANUP: Removed outdated patches for log2()/log10() and matrix().
#   They only applied to R v2.7.0 and before.
# 2010-02-10
# o Now also patches for R.filesets and R.utils are loaded, if available,
#   when aroma.core is loaded.
# 2009-09-04
# o Now smoothScatter() is copying the one in 'geneplotter' v1.2.4 or older,
#   if not R v2.9.0.
# 2009-05-13
# o Now the search() path is fixed for aroma.core as well.
# 2009-02-22
# o Now R.utils Settings object 'aromaSettings' is loaded/assign.
# 2008-07-24
# o Added patch for serialize() on Windows.
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
