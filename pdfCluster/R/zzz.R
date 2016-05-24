# .noGenerics <- TRUE
# .conflicts.OK <- TRUE

.onLoad <- function(lib, pkg)
{
    library.dynam("pdfCluster", pkg, lib)
   }

.onAttach <- function(lib, pkg)
{
    ver <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
    packageStartupMessage(paste(pkg, ver))
    packageStartupMessage( 
     "\n", 
       "PLEASE NOTE:  New features have been introduced in version 1.0-0", "\n",
       "These involve some changes in the package options", "\n",
       "see \"help(\"pdfCluster-package\")\" for the setting which", "\n",
       "reproduce the functioning of the previous versions.", " \n","\n")
}
