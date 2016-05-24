# init.R  
# Copyright (C) 2014 by Matthew Clegg

#  Initialization of the egcm module

#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

.onLoad <- function (libname, pkgname) {
#    packageStartupMessage("Simplified Engle-Granger Cointegration Models, version 1.0")
#    packageStartupMessage("Copyright (C) 2014 by Matthew Clegg.  Made available under GPL v2 or GPL v3")
    egcm.set.default.urtest("pp")
    egcm.set.default.i1test("pp")
    egcm.set.default.pvalue(0.05)
}