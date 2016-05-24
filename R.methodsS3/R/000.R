##############################################################################
# This code has to come first in a library. To do this make sure this file
# is named "000.R" (zeros).
##############################################################################

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# NAMESPACE: export()
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Sets attribute export to TRUE
export <- function(x) {
  attr(x, "export") <- TRUE;
  x;
}
export <- export(export)

# Sets attribute export to 'value'.
"export<-" <- export(function(x, value) {
  attr(x, "export") <- value;
  x;
})

noexport <- export(function(x) {
  attr(x, "export") <- FALSE;
  x;
})


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# NAMESPACE: S3method()
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Sets attribute 'S3class' to 'value'.
"S3class<-" <- export(function(x, value) {
  attr(x, "S3class") <- value;
  x;
})



# Use by setGenericS3() and setMethodS3()
.findFunction <- function(name, envir, inherits=rep(FALSE, times=length(envir))) {
  # Argument 'envir':
  if (!is.list(envir)) {
    envir <- list(envir);
  }
  n <- length(envir);

  # Argument 'inherits':
  inherits <- as.logical(inherits);
  stopifnot(length(inherits) == n);

  fcn <- pkg <- NULL;
  for (kk in seq_along(envir)) {
    env <- envir[[kk]];
    inh <- inherits[kk];
    if (exists(name, mode="function", envir=env, inherits=inh)) {
      fcn <- get(name, mode="function", envir=env, inherits=inh);
      pkg <- attr(env, "name");
      if (is.null(pkg)) {
        pkg <- "base"
        if (identical(env, baseenv())) {
        } else if (identical(env, globalenv())) {
          pkg <- "<R_GlobalEnv>"
        }
      } else {
        pkg <- gsub("^package:", "", pkg);
      }
      break;
    }
  } # for (kk ...)

  list(fcn=fcn, pkg=pkg);
} # .findFunction()


############################################################################
# HISTORY:
# 2013-10-06
# o Added .findFunction().
# 2012-04-17
# o Added S3class() function.
# o Added export() and noexport() functions.
# 2007-09-17
# o Removed support for R v2.2.0 and before by removing patch for missing
#   baseenv().
# 2007-04-07
# o Removed support for R v2.0.0 and before.
# 2006-02-09
# o Added baseenv() for R versions (< v2.2.0) where it does not exist.
#   This is used in setGenericS3() and setMethodS3() from R v2.3.0.
# 2005-02-15
# o Now require() is only called for R v1.9.1 or eariler.
# 2005-02-10
# o Moved R.KEYWORDS into its own source file.
# 2003-05-06
# o Added require(methods) to make sure getMethods() etc works.
# 2002-11-21
# o Added "..." to R.KEYWORDS.
# 2002-10-17
# o Removed obsolete "modifiers<-"().
# o Added also "Object" to the class attribute to make static methods to
#   work.
# 2002-10-16
# o There are times when
#     generic <- function(...) UseMethod()
#   is not working, for example
#     fcn <- get("generic"); fcn(myObj, ...);
#   For this reason, always do method dispatching using the name explicitly;
#     generic <- function(...) UseMethod("generic")
#
# 2002-10-15
# o Created from R.oo Object.R and ideas as described on
#    http://www.maths.lth.se/help/R/
############################################################################

