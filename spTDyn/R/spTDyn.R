####################################################################################


library(coda, warn.conflicts = FALSE)
library(sp, warn.conflicts = FALSE)
#library(forecast, warn.conflicts = FALSE)
#library(spacetime, warn.conflicts = FALSE)


.onLoad <-
    function(libname, pkgname)
{
      library.dynam(pkgname, pkgname, lib.loc=libname)
}


.onAttach <-
    function(libname, pkgname)
{
        ## figureout the version automatically
        library(help=spTDyn)$info[[1]] -> version
        version <- version[pmatch("Version",version)]
        um <- strsplit(version," ")[[1]]
        version <- um[nchar(um)>0][2]
        packageStartupMessage("\n## spTDyn version: ", version," \n")
}


####################################################################################
