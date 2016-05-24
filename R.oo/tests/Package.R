message("TESTING: Package...")

library("R.oo")

# By defining .onAttach() as follows in zzz.R for a package, an
# instance of class Package with the same name as the package will
# be made available on the search path. More over, the code below
# will also inform the user that the package has been loaded:
#
#  > library(R.oo)
#  R.oo v0.52 (2003/04/13) was successfully loaded.
#
.onAttach <- function(libname, pkgname) {
  pkg <- Package(pkgname);
  assign(pkgname, pkg, pos=getPosition(pkg));
  cat(getName(pkg), " v", getVersion(pkg), " (", getDate(pkg), ")",
    " was successfully loaded.\n", sep="");
}

# The Package class works for any packages, loaded or not.

# Some information about the base package
pkg <- Package("base")
print(pkg)
# [1] "Package: base v1.6.2 (NA) is loaded (pos=5). The official webpage
#      is NA and the maintainer is R Core Team <R-core@r-project.org>. The
#      package is installed in c:/PROGRA~1/R/rw1062/library/base/."
print(list.files(Package("base")$dataPath))

# Some information about the R.oo package
print(R.oo::R.oo)
# [1] "Package: R.oo v0.52 (2003/04/13) is loaded (pos=2). The official
#      webpage is http://www.braju.com/R/ and the maintainer is Henrik
#      Bengtsson <henrikb@braju.com>. The package is installed in
#      c:/PROGRA~1/R/rw1062/library/R.oo/."


pkg <- Package("R.oo")
classes <- getClasses(pkg)
print(classes)
stopifnot(all(c("Object", "Class", "Interface", "Exception", "Package") %in% classes))

pkg <- Package("R.oo")
res <- showDescriptionFile(pkg, pager=function(...) TRUE)
stopifnot(isTRUE(res))
res <- showNews(pkg, pager=function(...) TRUE)
stopifnot(isTRUE(res))
res <- showChangeLog(pkg, pager=function(...) TRUE)
stopifnot(isTRUE(res))
res <- showHistory(pkg, pager=function(...) TRUE)
stopifnot(isTRUE(res))


res <- try(Package("Non-Existing-Package"), silent=TRUE)
stopifnot(inherits(res, "try-error"))

message("TESTING: Package...DONE")
