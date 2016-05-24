
## Up until R 2.15.0, the require("methods") is needed but (now)
## triggers an warning from R CMD check
# .onLoad <- function(libname, pkgname){
#    #require("methods")  ## needed with R <= 2.15.0
#    loadRcppModules()
# }



## For R 2.15.1 and later this also works. Note that calling loadModule() triggers
## a load action, so this does not have to be placed in .onLoad() or evalqOnLoad().
Rcpp::loadModule("MODULE_DATA", TRUE)
Rcpp::loadModule("MODULE_PAR_KS", TRUE)


.onAttach <- function(libname, pkgname)
{
  packageStartupMessage("\n ClustMMDD = Clustering by Mixture models for Discrete Data.
  \n Version 1.0.3
  \n ClustMMDD is the R version of the stand alone c++ package named 'MixMoGenD'
  \n   that is available on www.u-psud.fr/math/~toussile.")
  
  packageStartupMessage("\n initializing ...", appendLF = FALSE)
  #Sys.sleep(1)
  packageStartupMessage(" Loaded \n")
}




