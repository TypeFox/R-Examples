##########################################################################
## start-up and clean-up functions
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## Copyright (C) 2011-present Ghislain Vieilledent
## 
##########################################################################

.onAttach <- function(...) {
   # echo output to screen
   packageStartupMessage("##\n## hSDM R package \n",
                         "## For hierarchical Bayesian species distribution models \n",
                         "##\n")
}

.onUnload <- function(libpath) {
    library.dynam.unload("hSDM", libpath)
}





