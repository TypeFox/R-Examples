'.onAttach' <- function(lib, pkg="gof")
  {    
    desc <- packageDescription(pkg)
    packageStartupMessage("Loading '", desc$Package, "' version ",desc$Version);
  }
