#
# 
# This file contains functions that deal with package lodaing and unloading.
# 

.onAttach <- function(libname, pkgname) {
  
#   # We load plyr with rollply as many functions will be directly useful to 
#   # interactive work with rollply (e.g. summarise)
#   if (!"package:plyr" %in% search()) {
#     packageStartupMessage('Loading required package: plyr')
#     library(plyr, quietly=TRUE)
#   }
# #   
#   # We load alphahull explicitely as it seems somehow badly designed: delvor
#   # does not find tri.pack() if the package is not explicitely loaded.
#   if (!"package:alphahull" %in% search()) {
#     packageStartupMessage('Loading required package: alphahull')
#     library(alphahull, quietly=TRUE)
#   }
  
}

















