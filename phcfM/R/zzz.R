##########################################################################
## start-up and clean-up functions
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 3, June 2007.  See the package LICENSE
## file for more information.
##
## Copyright (C) 2012-present Ghislain Vieilledent
## 
##########################################################################

.onAttach <- function(...) {
   # echo output to screen
   packageStartupMessage("##\n## phcfM R package \n",
                         "## A package for modelling anthropogenic deforestation \n",
                         "## Author: Ghislain Vieilledent <ghislain.vieilledent@cirad.fr> \n",
                         "## Support provided by the Cirad - http://www.cirad.fr/en \n##\n")
}

.onUnload <- function(libpath) {
    library.dynam.unload("phcfM", libpath)
}

