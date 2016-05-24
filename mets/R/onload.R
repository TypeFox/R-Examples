'.onAttach' <- function(libname, pkgname="mets")
  {    
    desc <- utils::packageDescription(pkgname)
    ## packageStartupMessage("Loading '", desc$Package, "' package...\n",
    ##                       "\tVersion: ", desc$Version, "\n",
    ##                       "\tOverview: help(package=", desc$Package, ")");
    packageStartupMessage(desc$Package, " version ",desc$Version);
  }
