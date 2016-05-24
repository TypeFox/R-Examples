## covr: skip=all

# Allows conflicts. For more information, see library() and
# conflicts() in [R] base.
.conflicts.OK <- TRUE;

.onLoad <- function(libname, pkgname) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Register vignette engines
  # NOTE: Are we doing this here because of backward compatibility
  #       issues, e.g. early version of R 3.0.x? /HB 2013-09-28
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  try({
    vignetteEngine <- get("vignetteEngine", envir=asNamespace("tools"));
    # Add vignette engine, iff missing.  But why?!? /HB 2013-09-28
    if (is.element("R.rsp::rsp", names(tools::vignetteEngine()))) {
      vignetteEngine("rsp", package=pkgname, pattern="[.][^.]*[.]rsp$",
                      weave=R.rsp::rspWeave, tangle=R.rsp::rspTangle);
    }
  }, silent=TRUE)
  ns <- getNamespace(pkgname);
  pkg <- Package(pkgname);
  assign(pkgname, pkg, envir=ns);

  # Set default global options when package is loaded
  devOptions("*", reset=TRUE)
}


.onAttach <- function(libname, pkgname) {
  startupMessage(get(pkgname, envir=getNamespace(pkgname)));
}


############################################################################
# HISTORY:
# 2013-05-30
# o Now registering RSP vignette engines.
# 2010-11-05
# o Created.
############################################################################
