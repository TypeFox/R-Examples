#.First.lib <- function(lib, pkg) library.dynam("sde", pkg, lib) 

#.noGenerics <- TRUE

.onAttach <- function(libname, pkgname)
{
    # require(methods)
 
    # require(zoo)
 packageStartupMessage(rep("#",44))
 packageStartupMessage("This is YUIMA Project package.")  
 packageStartupMessage("Check for the latest development version at:")
 packageStartupMessage("http://R-Forge.R-Project.org/projects/yuima")
 packageStartupMessage(rep("#",44))   
    
# require(KernSmooth, quietly=TRUE)
# library.dynam("yuima", pkgname, libname) 
}

