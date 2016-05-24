.onAttach <- function(libname, pkgname)
{
    packageStartupMessage("\n'drc' has been loaded.\n")
    packageStartupMessage("Please cite R and 'drc' if used for a publication,")
    packageStartupMessage("for references type 'citation()' and 'citation('drc')'.\n")
    
#    cat("\n")
#    cat("'drc' has been loaded.\n\n")
#    cat("for references type 'citation()' and 'citation('drc')'.\n\n")
}
