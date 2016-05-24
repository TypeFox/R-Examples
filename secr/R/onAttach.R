###############################################################################
## package 'secr'
## onAttach.R
## last changed 2011-06-16 2013-04-20
###############################################################################

.onAttach <- function (libname, pkgname) {
    version <- packageVersion('secr')
    packageStartupMessage( "This is secr ", version,
        ". For overview type ?secr" )
}
