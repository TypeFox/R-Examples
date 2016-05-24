
.onAttach <- function(lib, pkg) {
     packageStartupMessage(
         "This package is obsolete.\n",
         "The same functionality is in the R function 'psnice'\n",
         "    in the 'tools' package, which comes with R.\n")
}

