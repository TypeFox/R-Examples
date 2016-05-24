# Part of the mi package for multiple imputation of missing data
# Copyright (C) 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015 Trustees of Columbia University
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

.onLoad <- function(lib, pkg) {
  #   library.dynam("mi", pkg, lib)
  return(invisible(NULL))
}

.onUnload <- function(libpath) {
  #   library.dynam.unload("mi", libpath)
  return(invisible(NULL))
}

.onAttach <- function( ... ) {
  miLib <- dirname(system.file(package = "mi"))
  version <- utils::packageDescription("mi", lib.loc = miLib)$Version
  builddate <- utils::packageDescription("mi", lib.loc = miLib)$Packaged
  packageStartupMessage(paste("mi (Version ", version, ", packaged: ", builddate, ")", sep = ""))
  packageStartupMessage("mi  Copyright (C) 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015 Trustees of Columbia University")
  packageStartupMessage("This program comes with ABSOLUTELY NO WARRANTY.")
  packageStartupMessage("This is free software, and you are welcome to redistribute it")
  packageStartupMessage("under the General Public License version 2 or later.")
  packageStartupMessage("Execute RShowDoc('COPYING') for details.")
}
