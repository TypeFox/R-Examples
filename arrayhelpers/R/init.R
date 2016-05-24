.onAttach <- function (libname, pkgname){
  ## unlockBinding(".options", asNamespace("arrayhelpers")) 

  desc <- utils::packageDescription("arrayhelpers")
  vers <- paste("V. ", desc$Version)

  packageStartupMessage ("Package ",  desc$Package, ", version ", desc$Version, "\n\n",
       "If you use this package please cite it appropriately.\n",
       "   citation(\"", desc$Package, "\")\nwill give you the correct reference.", "\n\n",
       "The project homepage is ", desc$URL, "\n\n",
       sep = "")
}

