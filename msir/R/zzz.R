.onAttach <- function(lib, pkg)
{
  version <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
  packageStartupMessage("Package 'msir' version ", version)
  packageStartupMessage("Type 'citation(\"msir\")' for citing this R package in publications.")
  invisible()
}
  