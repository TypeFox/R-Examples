.onLoad <- function(libname, pkgname) {
    library.dynam("R2SWF", pkgname, libname);
    .Call("swfInit", PACKAGE = "R2SWF");
}

.onUnload <- function(libpath) {
    .C("Ming_collectGarbage", PACKAGE = "R2SWF");
    library.dynam.unload("R2SWF", libpath);
}

