## RShark -- R interface to the Shark libraries
##
## Copyright (C) 2010  Shane Conway	<shane.conway@gmail.com>
##
## This file is part of the RShark library for GNU R.
## It is made available under the terms of the GNU General Public
## License, version 2, or at your option, any later version,
## incorporated herein by reference.
##
## This program is distributed in the hope that it will be
## useful, but WITHOUT ANY WARRANTY; without even the implied
## warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
## PURPOSE.  See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public
## License along with this program; if not, write to the Free
## Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
## MA 02111-1307, USA


#' @useDynLib RcppShark
#' @importFrom Rcpp sourceCpp

 
# create own environment
RcppSharkEnv = new.env(parent = emptyenv())


 .onLoad <- function(libname, pkgname) {
  op <- options()
  op.devtools <- list(
    devtools.path = "~/R-dev",
    devtools.install.args = "",
    devtools.name = "Aydin Demircioglu",
    devtools.desc.author = '"Aydin Demircioglu <aydin.demircioglu@ini.rub.de> [aut, cre]"',
    devtools.desc.license = "GPL-3 + file LICENSE",
    devtools.desc.suggests = NULL,
    devtools.desc = list()
  )
  toset <- !(names(op.devtools) %in% names(op))
  if(any(toset)) options(op.devtools[toset])
  
#  find_rtools()
  
  invisible()
}


.onAttach <- function (libname, pkgname) {
        packageStartupMessage("RcppShark v0.1 loaded.")
}

