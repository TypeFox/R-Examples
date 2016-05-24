
## display version number and date when the package is loaded
#.onAttach <- function(libname, pkgname) {
#  desc  <- packageDescription(pkgname, libname)
#  packageStartupMessage(
#    'Package:  xergm.common\n', 
#    'Version:  ', desc$Version, '\n', 
#    'Date:     ', desc$Date, '\n', 
#    'Authors:  Philip Leifeld (Eawag and University of Bern)',
#    '\n\nPlease cite the xergm package in your publications ', 
#    '-- see citation("xergm").'
#  )
#}

# generics for extracting the formula from an estimation object
setGeneric("getformula", function(x) standardGeneric("getformula"), 
    package = "xergm.common")

# generic interpretation function
setGeneric("interpret", function(object, ...) standardGeneric("interpret"), 
    package = "xergm.common")

# generics for goodness-of-fit assessment
setGeneric("gof", function(object, ...) standardGeneric("gof"), 
    package = "xergm.common")
    
# generic checkdegeneracy function
setGeneric("checkdegeneracy", function(object, ...) 
    standardGeneric("checkdegeneracy"), package = "xergm.common")

