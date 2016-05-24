.onAttach <- function(libname, pkgname) {
    Pver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")
    packageStartupMessage(paste(pkgname, Pver))
    packageStartupMessage("Companion package to the book")
    packageStartupMessage(sQuote('Simulation and Inference for Stochastic Differential Equations With R Examples'))
    packageStartupMessage("Iacus, Springer NY, (2008)")
    packageStartupMessage("To check the errata corrige of the book, type vignette(\"sde.errata\")")
}

#.First.lib <- function(lib, pkg) library.dynam("sde", pkg, lib) 

#.noGenerics <- TRUE

#.onLoad <- function(libname, pkgname)
#{
# library.dynam("sde", pkgname, libname) 
#}
