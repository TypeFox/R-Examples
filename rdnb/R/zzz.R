.onAttach <- 
function(libname, pkgname) {
    ver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname), fields="Version")
    packageStartupMessage(" ")
    packageStartupMessage(paste("This is", pkgname, ver))
    packageStartupMessage(" ")
    packageStartupMessage("Type changes(\"rdnb\") to see changes/bug fixes, help(rdnb) for documentation")
    packageStartupMessage("or citation(\"rdnb\") for how to cite rdnb.")
    packageStartupMessage(" ")
}


#' @title View changes notes.
#' @description \code{changes} brings up the NEWS file of the package. 
#'
#' @param pkg Set to the default "rdnb". Other packages make no sense.
#' @export
#' @examples
#' \dontrun{
#' changes()
#' }
changes <- 
function(pkg="rdnb") {
    if(pkg=="rdnb") file.show(file.path(system.file(package="rdnb"), "NEWS"))
}
