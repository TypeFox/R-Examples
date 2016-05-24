.onLoad <- function(lib, pkg) {
library.dynam("seewave", pkg, lib)
}

## .onAttach <- function(...){
## packageStartupMessage(paste("
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Welcome to seewave ! 
## The package is regularly updated, please check for new version [http://rug.mnhn.fr/seewave]
## Thanks to use the right reference when citing seewave in publications
## See citation('seewave')
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## "))
## }

