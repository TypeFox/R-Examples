##
## R source file
## This file is part of rgl
##
## $Id: plugin.R 376 2005-08-03 23:58:47Z dadler $
##

##
## quit R plugin
## 
##

rgl.quit <- function() {

  unloadNamespace("rgl")

}


