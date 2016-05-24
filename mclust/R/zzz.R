.onLoad <- function(libname, pkgname) 
{
  library.dynam("mclust", pkgname, libname)
}

.onAttach <- function(lib, pkg)
{
  unlockBinding(".mclust", asNamespace("mclust")) 
  version <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
  
  if(interactive())
    { # > figlet -f slant MCLUST
      packageStartupMessage(
"    __  ___________    __  _____________
   /  |/  / ____/ /   / / / / ___/_  __/
  / /|_/ / /   / /   / / / /\\__ \\ / /   
 / /  / / /___/ /___/ /_/ /___/ // /    
/_/  /_/\\____/_____/\\____//____//_/    version ", version)
}
else
  { packageStartupMessage("Package 'mclust' version ", version) } 

  packageStartupMessage("Type 'citation(\"mclust\")' for citing this R package in publications.")
  invisible()
}



