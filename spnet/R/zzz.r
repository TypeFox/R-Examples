#.onLoad <- function(libname, pkgname) {
#    if (!require(methods)) {
#        stop("Require methods package")
#    }
#}



.onAttach <- function(libname, pkgname) {
  packageStartupMessage('\n', 'Welcome to spnet', '.')
  packageStartupMessage('You are running version ', utils::packageVersion("spnet"), '.\n')
  packageStartupMessage(
    'For introductory material, type ',
    #"'vignette(package=\"spnet\")'.\n"
    "`vignette('spnet-overview')`.\n"
  )
  packageStartupMessage(
    'If you use this package in your work, thank you for rewarding our work by citing the package in your own one. ',
    "Please type `citation(\'spnet\')` for citation information.\n"
  )
  packageStartupMessage(
    "Type `?spnet` ",
    "to display the help index.\n"
  )
  
}

#' @export
.Last.lib <- function(libpath) {
  message('\n', "Thank you for using the 'spnet' package", '.')
  message('See you soon', '!')
}

.onUnload <- function(libpath) {
  #library.dynam.unload("spnet", libpath )
}
