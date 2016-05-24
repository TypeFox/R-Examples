.onLoad <- function(libname, pkgname)
{ 
    # Load dll:
    library.dynam("CAMAN", pkgname, libname)
}