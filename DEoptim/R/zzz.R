".onLoad" <- function (lib, pkg)
{
  library.dynam(pkg, pkg, lib)
  
}
".onAttach" <- function (lib, pkg)
{
  packageStartupMessage("\nDEoptim package",
   "\nDifferential Evolution algorithm in R",
   "\nAuthors: D. Ardia, K. Mullen, B. Peterson and J. Ulrich\n")
}
