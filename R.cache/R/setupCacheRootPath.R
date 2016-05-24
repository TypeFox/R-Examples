#########################################################################/**
# @RdocDefault setupCacheRootPath
#
# @title "Interactively offers the user to set up the default root path"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{defaultPath}{Default root path to set.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns (invisibly) the root path,
#   or @NULL if running a non-interactive session.
# }
#
# \details{
#   If the cache root path is already set, it is used and nothing is done.
#   If the "default" root path (\code{defaultPath}) exists, it is used,
#   otherwise, if running interactively, the user is asked to approve
#   the usage (and creation) of the default root path.
#   In all other cases, the cache root path is set to a session-specific
#   temporary directory.
# }
#
# @author
#
# \seealso{
#  Internally, @see "setCacheRootPath" is used to set the cache root path.
#  The @see "base::interactive" function is used to test whether \R is
#  running interactively or not.
# }
#
# @keyword "programming"
# @keyword "IO"
# @keyword "internal"
#*/#########################################################################
setMethodS3("setupCacheRootPath", "default", function(defaultPath="~/.Rcache/", ...) {
  rootPath <- getCacheRootPath(NULL);

  # If already set, nothing to do.
  if (!is.null(rootPath)) {
    return(invisible(rootPath));
  }

  # Use a temporary root path...
  rootPath <- file.path(tempdir(), ".Rcache");

  # unless the default directory exists, ...
  if (isDirectory(defaultPath)) {
    rootPath <- defaultPath;
  } else if (interactive()) {
    # or we cn ask the user to confirm the default path...
    prompt <- "The R.cache package needs to create a directory that will hold cache files.";
    if (identical(defaultPath, "~/.Rcache/")) {
      prompt <- c(prompt, "It is convenient to use one in the user's home directory, because it remains also after restarting R.");
    }
    prompt <- c(prompt, sprintf("Do you wish to create the '%s' directory? If not, a temporary directory (%s) that is specific to this R session will be used.", defaultPath, rootPath));
    prompt <- paste(prompt, collapse=" ");
    tryCatch({
      ans <- .textPrompt(prompt=prompt, options=c("Y"="yes", "n"="no"));
      if (ans == "yes") rootPath <- defaultPath;
    }, condition=function(ex) {});
  }

  setCacheRootPath(rootPath);
  rootPath <- getCacheRootPath();

  invisible(rootPath);
}) # setupCacheRootPath()


############################################################################
# HISTORY:
# 2011-12-29
# o Added setupCacheRootPath().
############################################################################
