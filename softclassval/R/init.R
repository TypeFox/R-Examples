.onAttach <- function (libname, pkgname){
  ## unlockBinding(".options", asNamespace("softclassval")) 

  desc <- utils::packageDescription("softclassval")
  vers <- paste("V. ", desc$Version)

  packageStartupMessage ("Package ",  desc$Package, ", version ", desc$Version, "\n\n",
'Use of this package is granted under GPL Version 3 or above.
   RShowDoc ("GPL-3")
will display the license.

In addition, the use of ', desc$Package ,' must be cited appropriately.
   citation("', desc$Package, '")
will give you the correct reference.

The project homepage is ', desc$URL, '.

',
       sep = "")
}

