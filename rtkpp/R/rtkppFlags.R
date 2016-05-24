#-----------------------------------------------------------------------
#     Copyright (C) 2014  Serge Iovleff, University Lille 1, Inria
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
# Copyright (C) 2009 - 2013 Dirk Eddelbuettel and Romain Francois
#
# This file is part of Rcpp.
#
# Rcpp is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# Rcpp is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Rcpp.  If not, see <http://www.gnu.org/licenses/>.
#-----------------------------------------------------------------------
#' CxxFlags defaults for the rtkpp
#' @rdname rtkppFlags
#' @keywords internal
CxxFlags <- function(cpp11=FALSE) { cat(.rtkppCxxFlags(cpp11=cpp11)) }
#-----------------------------------------------------------------------
#' CppFlags defaults for the rtkpp
#' @rdname rtkppFlags
#' @keywords internal
CppFlags <- function() { cat(.rtkppCppFlags()) }
#' LdFlags defaults for the rtkpp package
#' @rdname rtkppFlags
#' @keywords internal
LdFlags <- function() { cat(.rtkppLdFlags()) }

# make sure system.file returns an absolute path
###########################
# Adapted from Rcpp package
###########################
# @rdname rtkppFlags
# @keywords internal
.rtkpp.system.file <- function(...)
{ tools::file_path_as_absolute( base::system.file( ..., package = "rtkpp" ) )}

# Provide compiler flags -- i.e. -I/path/to/RTKpp.h
# @keywords internal
.rtkppCxxFlags <- function(cpp11=FALSE)
{
  path1 <- .rtkpp.system.file( "include" )
  path2 <- .rtkpp.system.file( "projects" )
  if (.Platform$OS.type=="windows")
  {
    path1 <- .asBuildPath(path1);
    path2 <- .asBuildPath(path2);
  }
  paste("-I", path1, " -I", path2, if (cpp11) " -std=c++11 " else "", sep="")
}

# Provide internal STK++ macros
# @keywords internal
.rtkppCppFlags <- function(cpp11=FALSE)
{
  paste("-DIS_RTKPP_LIB -DSTKUSELAPACK", sep="")
}

# Provide linker flags -- i.e. /path/to/rtkpp.so
###########################
# Adapted from Rcpp package
###########################
# @rdname rtkppFlags
# @keywords internal
.rtkppLdFlags <- function()
{
  if (nzchar(.Platform$r_arch))
  {	## eg amd64, ia64, mips
    path <- .rtkpp.system.file("lib",.Platform$r_arch)
  }
  else { path <- .rtkpp.system.file("lib") }
  paste( path, "/libSTKpp.a", sep="")
}

# Transform a path for passing to the build system on the command line.
# Leave paths alone for posix. For Windows, mirror the behavior of the
# R package build system by starting with the fully resolved absolute path,
# transforming it to a short path name if it contains spaces, and then
# converting backslashes to forward slashes
###########################
# Adapted from Rcpp package
###########################
# @rdname rtkppFlags
# @keywords internal
.asBuildPath <- function(path)
{
  if (.Platform$OS.type == "windows")
  {
    path <- normalizePath(path)
    if (grepl(' ', path, fixed=TRUE))
      path <- utils::shortPathName(path)
    path <- gsub("\\\\", "/", path)
  }
  return(path)
}
