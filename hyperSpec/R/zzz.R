.onAttach <- function (libname, pkgname){
  unlockBinding(".options", asNamespace("hyperSpec"))

  desc <- utils::packageDescription("hyperSpec")
  vers <- paste("V. ", desc$Version)

  packageStartupMessage ("Package ",  desc$Package, ", version ", desc$Version, "\n\n",
       "To get started, try\n",
       '   vignette ("introduction", package = "hyperSpec")\n',
       '   package?hyperSpec \n',
       '   vignette (package = "hyperSpec")\n\n',
       "If you use this package please cite it appropriately.\n",
       "   citation(\"hyperSpec\")\nwill give you the correct reference.", "\n\n",
       "The project homepage is http://hyperspec.r-forge.r-project.org\n\n",
       sep = "")
}

