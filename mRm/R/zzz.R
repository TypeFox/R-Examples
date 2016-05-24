#    Copyright (C) 2012  David Preinerstorfer
#    david.preinerstorfer@univie.ac.at
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details. A copy may be obtained at
#    http://www.r-project.org/Licenses/

.onAttach <- function(...) {

  packageStartupMessage("## mRm - version 1.1.5 - License:", " ", "GPL-2 \n", sep="")
  packageStartupMessage("## NO WARRANTY PROVIDED")

}

.onUnload <- function(libpath) {
    library.dynam.unload("mRm", libpath)
}

