## RcppShark.package.skeleton.R: makes a skeleton for a package that wants to use RcppShark
##
## Copyright (C)  2014  Qiang Kou
##
## This file was part of RcppMLPACK.
##
## RcppMLPACK is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## RcppMLPACK distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with RcppMLPACK.  If not, see <http://www.gnu.org/licenses/>.

inlineCxxPlugin <- function(...) {
    plugin <- Rcpp::Rcpp.plugin.maker(
#        include.before = "#include <sharkHeader.h>", 
        libs = sprintf("%s $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)", RcppSharkLdFlags()),
        LinkingTo      = c( "Rcpp", "BH", "RcppShark"),
        package        = "RcppShark"
    )
    settings <- plugin()
    settings
}
