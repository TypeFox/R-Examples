.noGenerics <- TRUE

.onLoad <- function(libname, pkgname) 
{
    chkernel()
    library.dynam("spatialkernel", pkgname, libname)
}

.onAttach <- function(libname, pkgname)
{
    desfile <-  file.path(libname, pkgname, "DESCRIPTION",
                          fsep = .Platform$file.sep)
    verline <- readLines(desfile, n=2, ok=TRUE)[[2]]
    version <- gsub("Version: ", "", verline, ignore.case=TRUE)
    packageStartupMessage("\nThis is ", pkgname," ", version, "\n\n")
}

.onUnload <- function(libpath)
    library.dynam.unload("spatialkernel", libpath)
