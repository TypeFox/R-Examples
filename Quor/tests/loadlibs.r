## File Quor/tests/loadlibs.r
##
## Quor package for R (http://www.R-project.org)
## Copyright (C) 2014 Adriano Polpo, Carlos A. de B. Pereira, Cassio P. de Campos.
##
##    This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program.  If not, see <http://www.gnu.org/licenses/>.

## This code is used for testing purposes. The Quor library does not
## depend on it for any of its functionalities

dependencies <- function() {
    return(c("combinat"))
}

installpacks <- function(loc=NULL,repos="http://stat.ethz.ch/CRAN/") {
  ## set the repository to use
  options(repos=repos)
  ## install the packages
  
  for(pack in dependencies()) {
      install.packages(pack,lib=loc)
  }
  
  ## this following line install the Quor package itself, so nothing else is needed.
  ## For testing, sometimes it is better to work without installing it for a while...
  ##      install.packages('./Quor_version.tar.gz',repos=NULL,type="source")
}


loadlibs <- function(libdir=NULL,package.loc='..',forcelocal=FALSE) {
  w <- options("warn")
  options("warn" = -1)
  if (forcelocal==TRUE || require("Quor",quietly=TRUE)==FALSE) {
    for(pack in dependencies()) {
        library(package=pack,lib.loc=libdir,character.only=TRUE)
    }
    for(fn in dir(file.path(package.loc,"R"),pattern="[.]r$")) {
        source(file.path(package.loc,"R",fn))
    }
    if(is.loaded('quorccore')) {
        dyn.unload(file.path(package.loc,'src','Quor.so'))
    }
    system(paste("R CMD SHLIB -o ",file.path(package.loc,"src","Quor.so")," ",file.path(package.loc,"src","core.c"),sep=""))
    dyn.load(file.path(package.loc,'src','Quor.so'))
  } else {
    library("Quor")
  }
  options("warn" = w[[1]])
}
