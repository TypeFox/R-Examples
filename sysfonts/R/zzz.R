.onLoad <- function(libname, pkgname) {
    library.dynam("sysfonts", pkgname, libname);
    .add.default.font.paths();
    .add.default.fonts();
}

.onUnload <- function(libpath) {
    .clean.fonts();
    library.dynam.unload("sysfonts", libpath);
}

