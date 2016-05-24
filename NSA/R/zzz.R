# Allows conflicts. For more information, see library() and
# conflicts() in [R] base.
.conflicts.OK <- TRUE;

# WORKAROUND: In order for the package to work with the most recent
# version of R devel, which automatically add namespaces to packages
# who do not have one, we explicitly have specify the following.
# /HB 2011-08-01
# From R.utils:
cat <- R.utils::cat;


.First.lib <- function(libname, pkgname) {
  pd <- packageDescription(pkgname);

  packageStartupMessage(pkgname, " v", pd$Version, " (", 
    pd$Date, ") successfully loaded. See ?", pkgname, " for help.");
}
