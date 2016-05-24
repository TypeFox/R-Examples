################################################################################
######################################### zzz.R ################################
################################################################################

## Display a startup message

.onAttach <- function(libname, pkgname) {

  packageStartupMessage('This is relaxnet version ',
                        packageDescription('relaxnet')$Version, "\n")

}
