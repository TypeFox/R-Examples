#-----------------------------------------------------------------------
#     Copyright (C) 2012-2014  Serge Iovleff, University Lille 1, Inria
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as
#    published by the Free Software Foundation; either version 2 of the
#    License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public
#    License along with this program; if not, write to the
#    Free Software Foundation, Inc.,
#    59 Temple Place,
#    Suite 330,
#    Boston, MA 02111-1307
#    USA
#
#    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
#
#-----------------------------------------------------------------------
# plugins are functions that return a list with :
#
# includes
# mandatory. it is included at the top of the compiled file by cxxfunction
#
# body
# optional. a function that takes one argument (the body of the c++ function)
# and returned a modified version of the body. The "Rcpp" plugin uses this to
# surround the code with the BEGIN_RCPP and END_RCPP macros
#
# LinkingTo
# optional. character vector containing the list of packages that the code
# needs to link to. This adds the include path of the given packages. The "Rcpp"
# and "RcppArmadillo" plugins use this.
#
# env
# optional. named list of environment variables. For example, the "Rcpp"
# plugin uses this to add Rcpp user library to the PKG_LIBS environment variable.
#
#
# plugins can be manually registered using the registerPlugin function.
# Alternatively, a package may supply an inline plugin implicitely by defining a
# function called inlineCxxPlugin, which does not necessarily need to be
# exported from the namespace of the package.
inlineCxxPlugin <- function()
{
includes <- "
#include <RTKpp.h>
#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

using namespace Rcpp;
using namespace STK;
"

body <- function( x ){ sprintf( "BEGIN_RCPP\n%s\nEND_RCPP", x ) }

LinkingTo <- unique(c("Rcpp", "rtkpp" ))

#  libs <-""
libs <- "`$(R_HOME)/bin/Rscript -e \"rtkpp:::LdFlags()\"` `$(R_HOME)/bin/Rscript -e \"Rcpp:::LdFlags()\"`"
Cxx = "`${R_HOME}/bin/Rscript -e \"rtkpp:::CxxFlags()\"` `$(R_HOME)/bin/Rscript -e \"Rcpp:::CxxFlags()\"` -DIS_RTKPP_LIB -DSTKUSELAPACK"

out <- list( env = list( PKG_LIBS = libs, PKG_CXXFLAGS = Cxx ), includes = includes, LinkingTo = LinkingTo , body = body)

out
}

