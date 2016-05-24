################################################################################
######################################### zzz.R ################################
################################################################################

## Display a startup message

.onAttach <- function(libname, pkgname) {

  packageStartupMessage('This is widenet version ',
                        packageDescription('widenet')$Version, "\n")

}
