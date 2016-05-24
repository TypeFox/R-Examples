## covr: skip=all

# Detach the 'R.oo' attached in file 030.ObjectClassFunctions.R
if (is.element("R.oo", search())) detach("R.oo");

.onUnload <- function(libpath) {
##  message("R.oo::.onUnload()");
  # Force finalize() on Object:s
  base::gc();
} # .onUnload()


.onLoad <- function(libname, pkgname) {
  ns <- getNamespace(pkgname);

  ## Doing assign(pkgname, Package(pkgname), envir=ns) seems to
  ## introduce potential cyclic loading of the R.oo namespace.
  ## My best guess is that it has to do with garbage collection.
  ## Because of this, we use a "delayed" assignment. /HB 2013-10-10
  delayedAssign(pkgname, Package("R.oo"), eval.env=ns, assign.env=ns);

  # Create getCall() generic function, iff missing (R < 2.14.0)
  if (!exists("getCall", mode="function")) {
    assign("getCall", function(...) UseMethod("getCall"), envir=ns);
  }
} # .onLoad()


.onAttach <- function(libname, pkgname) {
  pkg <- get(pkgname, envir=getNamespace(pkgname));
  startupMessage(pkg);
} # .onAttach()



############################################################################
# HISTORY:
# 2014-02-21
# o Added .onUnload() which calls the garbage collector.
############################################################################
