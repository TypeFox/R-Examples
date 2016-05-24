'.onAttach' <- function(lib, pkg="gap")
  {   
    desc <- packageDescription(pkg)
    ## packageStartupMessage("Loading '", desc$Package, "' package...\n",
    ##                       "\tVersion: ", desc$Version, "\n",
    ##                       "\tOverview: help(package=", desc$Package, ")");
    packageStartupMessage(desc$Package, " version ",desc$Version);
    
  }

.noGenerics <- TRUE
.onUnload <- function(libpath) library.dynam.unload("gap", libpath)
