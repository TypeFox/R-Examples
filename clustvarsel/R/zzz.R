.onAttach <- function(lib, pkg)
{
  .mclust <- get(".mclust", envir = asNamespace("mclust"))
  .mclust$hcUse <- "SVD"
  assign(".mclust", .mclust, envir = asNamespace("mclust"))  
  
  version <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
  packageStartupMessage("Package 'clustvarsel' version ", version)
  packageStartupMessage("Type 'citation(\"clustvarsel\")' for citing this R package in publications.")
  invisible()
}



  