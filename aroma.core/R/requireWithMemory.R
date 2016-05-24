############################################################################
# requireWithMemory() works like require(package).  However, with for
# instance package="GLAD" the latter will open a modal GUI dialog if there
# are problems with the shared library which each time has to be closed by
# the user [1].  Since for instance the plotting of ChromosomeExplorer will
# keep trying to load GLAD that becomes super annoying.  This function will
# remember the failure and not try to load GLAD until 24 hours later,
# or if GLAD is reinstalled.
#
# REFERENCES:
# [1] Bioconductor developers list, thread 'GLAD: Suggestion: Throw a formal
#     error condition instead of modal GUI dialog', started on 2010-11-23.
############################################################################
setMethodS3("requireWithMemory", "default", function(package="GLAD", ..., since=24*3600, force=FALSE) {
  pkg <- "R.utils";
  require(pkg, character.only=TRUE) || throw("Package not loaded: ", pkg);

  # If package is already loaded, do nothing
  if (isPackageLoaded(package)) {
    return(TRUE);
  }

  # Argument 'package'
  package <- Arguments$getCharacter(package);

  # Argument 'since'
  since <- Arguments$getInteger(since, range=c(0, Inf));

  # Argument 'force'
  force <- Arguments$getLogical(force);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check cache for loading errors within the last 'since' seconds
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pd <- packageDescription(package);
  key <- list(package=package, packageDescription=pd);
  dirs <- c("aroma.core", "requireWithMemory", package);
  if (!force) {
    lastResult <- loadCache(key=key, dirs=dirs);
    if (!is.null(lastResult)) {
      # Too old? Try only once every 'since' seconds
      if (Sys.time() - lastResult$timestamp < since) {
        return(FALSE);
      }
    }
  }

  lastError <- NULL;
  errorMessage <- NULL;
  con <- textConnection("errorMessage", open="w", local=TRUE);
  sink(con, type="message");
  on.exit({
   if (!is.null(con)) {
     sink(type="message");
     close(con);
   }
  });
  tryCatch({
    library(package, character.only=TRUE, logical.return=FALSE);
  }, error = function(ex) {
    lastError <<- ex;
  });

  sink(type="message");
  close(con);
  con <- NULL;

  # Was the package loaded?
  res <- isPackageLoaded(package);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # If failing to load shared library, then memoize the failure
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # If not, then check what should happend next time
  if (!res) {
    # If cannot load shared library, then don't try again
    # for a while.
    errorMessage <- paste(errorMessage, collapse="\n");
    if (regexpr("LoadLibrary failure", errorMessage) != -1 ||
        regexpr("unable to load shared object", errorMessage) != -1) {
      lastResult <- list(
        timestamp = Sys.time(),
        errorMessage = errorMessage,
        error = lastError
      );
      # Save to cache
      saveCache(lastResult, key=key, dirs=dirs);
    }
  }

  res;
}, private=TRUE) # requireWithMemory()


############################################################################
# HISTORY:
# 2010-12-07
# o Added argument 'force' to requireWithMemory().
# 2010-12-06
# o Added requireWithMemory().
# o Created.
############################################################################
