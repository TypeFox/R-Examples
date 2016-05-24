###########################################################################/**
# @RdocFunction .makeObjectFinalizer
#
# @title "Creates a standalone finalizer function for Object"
#
# \description{
#   @get "title", which assumes only that the \pkg{base} package is available.
# }
#
# @synopsis
#
# \arguments{
#   \item{this}{The @Object to be finalized.}
#   \item{reloadRoo}{If @TRUE, the finalizer will try to temporary
#     reload the \pkg{R.oo} package if not loaded.}
# }
#
# \value{
#   Returns a @function that can be passed to @see "base::reg.finalizer".
# }
#
# \details{
#   The created finalizer is reentrant.
#   This is always the case when the \pkg{R.oo} package is already
#   loaded when the finalizer is called.
#   It is also always the case on R v2.15.2 Patched r61487 and beyond.
#   In other cases, the finalizer inspects the call stack
#   (via @see "sys.calls") to check whether @see "base::parse"
#   has been called or not.  If it is on the call stack, it indicates
#   that @see "base::parse" triggered the garbage collection, and
#   the \pkg{R.oo} package will \emph{not} be reloaded in order to
#   avoid risking @see "base::parse" being called again.
# }
#
# @keyword internal
#*/###########################################################################
.makeObjectFinalizer <- function(this, reloadRoo=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  getObjectInfo <- function(this) {
    env <- attr(this, ".env");
    if (is.null(env)) return(NA_character_);

    # base::environmentName() was added to R v2.5.0
    if (exists("environmentName", mode="function", envir=getNamespace("base"), inherits=FALSE)) {
      name <- environmentName(env);
    } else {
      name <- "";
    }

    if (name == "") {
      # Cannot assume 'utils' is loaded
      name <- utils::capture.output(print.default(env));
      name <- name[1L]; # Just in case
      name <- gsub("[<]*environment:[ ]*([^>]*)[>]", "\\1", name);
    }

    className <- class(this)[1L];
    name <- paste(className, ": ", name, sep="");
    # More debug information
    if (className == "Package") {
      name <- sprintf("%s {.name='%s'}", name, env$.name);
    }

    name;
  } # getObjectInfo()

  isLibraryReentrant <- function() {
    getRversion2 <- function() {
      rVer <- R.version[c("major", "minor", "status", "svn rev")];
      names(rVer)[3:4] <- c("patched", "rev");
      rVer$patched <- ifelse(identical(rVer$patched, "Patched"), 1L, 0L);
      rVer$rev <- ifelse(is.null(rVer$rev), 0L, rVer$rev);
      rVer <- lapply(rVer, FUN=as.numeric);
      rVer;
    } # getRversion2()

    rVer <- getRversion2();
    if (rVer$major >= 3) return(TRUE);
    if (rVer$major < 2) return(FALSE);
    if (rVer$minor >= 15.3) return(TRUE);
    if (rVer$minor < 15.2) return(FALSE);
    if (rVer$patched < 1) return(FALSE);
    if (rVer$rev < 61487) return(FALSE);
    TRUE;
  } # isLibraryReentrant()

  isParseCalled <- function() {
    calls <- sys.calls();
    if (length(calls) == 0L) return(FALSE);
    for (kk in seq_along(calls)) {
      call <- calls[[kk]];
      name <- call[[1L]];
      if (!is.symbol(name)) next;
      if (name == "do.call") {
        name <- call[[2L]];
        if (!is.symbol(name)) next;
        name <- as.character(name);
      }
      if (name == "parse") {
        return(TRUE);
      }
    } # for (kk ...)
    FALSE;
  } # isParseCalled()

  isLibraryActive <- function() {
    # Just in case the below won't work one day due to R updates...
    tryCatch({
      # Identify the environment/frame of interest by making sure
      # it at least contains all the arguments of source().
      argsToFind <- names(formals(base::library));

      # Scan the call frames/environments backwards...
      srcfileList <- list();
      for (ff in sys.nframe():0) {
        env <- sys.frame(ff);

        # Does the environment look like a library() environment?
        exist <- sapply(argsToFind, FUN=exists, envir=env, inherits=FALSE);
        if (!all(exist)) {
          # Nope, then skip to the next one
          next;
        }

        return(TRUE);
      } # for (ff ...)
    }, error = function() {});

    FALSE;
  } # isLibraryActive()

  isRooLoading <- function() {
    for (n in sys.nframe():1) {
      env <- sys.frame(n);
      if (exists("__NameSpacesLoading__", envir=env, inherits=FALSE)) {
        pkgs <- get("__NameSpacesLoading__", envir=env, inherits=FALSE);
        if (is.element("R.oo", pkgs)) return(TRUE);
      }
    } # for (n ...)
    FALSE;
  } # isRooLoading()

  isRooLoaded <- function() {
    is.element("R.oo", loadedNamespaces());
  } # isRooLoaded()

  # NOTE: The finalizer() depends on the 'this' object. # /HB 2011-04-02
  # Note, R.oo might be detached when this is called!  If so, reload
  # it, this will be our best chance to run the correct finalizer(),
  # which might be in a subclass of a different package that is still
  # loaded.
  finalizer <- function(env) {
    debug <- FALSE;
    if (debug) message(sprintf("finalizer(): reloadRoo=%s", reloadRoo));
    if (debug) message(paste(capture.output(print(sessionInfo())), collapse="\n"))
    # Classes for which it is known that finalizer() does nothing
    if (is.element(class(this)[1L], c("Object", "Class", "Package"))) {
      return();
    }

    # Do nothing if R.oo is being loaded
    if (isRooLoading()) return();

    # Is R.oo::finalize() available?
    if (isRooLoaded()) {
      if (debug) message(sprintf("finalizer(): Call finalize()..."));
      R.oo::finalize(this);
      if (debug) message(sprintf("finalizer(): Call finalize()...done"));
      return();
    }

    # Nothing to do?
    if (!reloadRoo) return();

    warning(sprintf("Object was not be finalize():d properly because the R.oo package was not loaded: %s", getObjectInfo(this)));
    return();

    # Skip the finalizer()?
    genv <- globalenv();
    if (exists("...R.oo::skipFinalizer", envir=genv, inherits=FALSE)) {
      if (debug) message(sprintf("Skipping finalizer()."));
      return();
    }

    # Skip because library() is active?
    if (isLibraryActive()) {
      warning(sprintf("Detected finalization of Object while library() is working on attaching a package and the R.oo package is not loaded.  This is most likely due to a *temporary* Object allocated in .onLoad()/.onAttach(), which is not allowed. Finalization of Object will *not* be skipped: %s", getObjectInfo(this)));
      return();
    }

    # Check if base::library() is reentrant...
    if (!isLibraryReentrant()) {
      # If not, check if base::parse() triggered the garbage collection
      # and/or has been called, because then we must not call library(),
      # because it will in turn call parse() potentially causing R to
      # crash.
      if (isParseCalled()) {
        warning("Object may not be finalize():d properly because the R.oo package was not loaded and will not be reloaded, because if done it may crash R (running version of R is prior to R v2.15.2 Patched r61487 and the garbage collection was triggered by base::parse()): ", getObjectInfo(this));
        return();
      }
    }

    warnMsg <- "";
    # (1) Load the 'R.oo' package
    if (debug) message(sprintf("finalizer(): Reloaded R.oo..."));
    suppressMessages({
      warnMsg <- suppressWarnings({
        loadNamespace("R.oo");
      });
    });
    if (debug) message(sprintf("finalizer(): Reloading R.oo...done"));
    if (is.character(warnMsg)) {
      warnMsg <- sprintf(" (with warning %s)", warnMsg);
    } else {
      warnMsg <- "";
    }

    # For unknown reasons R.oo might not have been loaded.
    if (!isRooLoaded()) {
      warning(sprintf("Object may not be finalize():d properly because the R.oo package failed to reload%s: %s", warnMsg, getObjectInfo(this)));
      return();
    }

    # Make sure to unload R.oo at the end
    on.exit(unloadNamespace("R.oo"));

    if (debug) message(sprintf("finalizer(): Call finalize()..."));
    R.oo::finalize(this);
    if (debug) message(sprintf("finalizer(): Call finalize()...done"));

    # NOTE! Before unloading R.oo again, we have to make sure the Object:s
    # allocated by R.oo itself (e.g. an Package object), will not reload
    # R.oo again when being garbage collected, resulting in an endless
    # loop.  We do this by creating a dummy finalize() function, detach
    # R.oo, call garbage collect to clean out all R.oo's objects, and
    # then remove the dummy finalize() function.
    # (2) Put a dummy finalize() function on the search path.
    if (!exists("...R.oo::skipFinalizer", envir=genv, inherits=FALSE)) {
      assign("...R.oo::skipFinalizer", TRUE, envir=genv, inherits=FALSE);
      on.exit({
        rm(list="...R.oo::skipFinalizer", envir=genv, inherits=FALSE);
      }, add=TRUE)
    }

    # (4) Force all R.oo's Object:s to be finalize():ed.
    if (debug) message(sprintf("finalizer(): Calling base::gc()"));
    base::gc();

    if (debug) message(sprintf("finalizer(): done"));
  } # finalizer()

  return(finalizer);
} # .makeObjectFinalizer()


############################################################################
# HISTORY:
# 2014-02-21
# o Now .makeObjectFinalizer() returns a finalizer that only loads the
#   R.oo package (it used to attach it).
# o Now calling base::gc() explicitly.
# 2013-10-13
# o ROBUSTNESS: Now Object finalizers will only try to re-attach the
#   'R.oo' package if library() is currently attaching a package.  This
#   can occur if a temporary Object is allocated in .onAttach().
# o Now warnings generated by the finalizer function returned by
#   .makeObjectFinalizer() is more informative on Package objects.
# 2013-09-20
# o BUG FIX: The finalizer returned by .makeObjectFinalizer() assumed
#   that the 'utils' is attached while calling capture.output(), which
#   under certain conditions could generate 'Error in getObjectInfo(this) :
#   could not find function "capture.output"'.
# 2013-01-08
# o ROBUSTNESS: Now .makeObjectFinalizer() returns a finalizer that is
#   reentrant, i.e. it will only try to reload R.oo on R versions where
#   library() is reentrant or when the garbage collector was not triggered
#   by base::parse(), otherwise it will not finalize the Object.
# o CLEANUP: Added internal .makeObjectFinalizer().
############################################################################
