if (version$language == "R") {
".onLoad" <- 
function(libname, pkgname) {
        .jpackage(pkgname, lib.loc=libname)
}
} else {
".onLoad" <- 
function(libname, pkgname) {
        .JavaAttachClassPath(system.file("java", "sjdbc.jar", package=pkgname))
}
}
