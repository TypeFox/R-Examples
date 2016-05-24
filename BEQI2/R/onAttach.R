.onAttach <-
function(libname, pkgname) {
    packageStartupMessage(
        "For more information on the usage of the BEQI2 tool, type: ", 
        'vignette("BEQI2")'
    )
}  
