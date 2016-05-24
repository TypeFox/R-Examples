## ----------------------------------------------------------------------------
## This file is part of boolean3
##
## Copyright (C) 2011--2014 Jason W. Morgan <morgan.746@osu.edu>
##
## boolean3 represents a substantial re-write of the original boolean package
## developed by Bear Braumoeller, Ben Goodrich, and Jacob Kline. This version
## was developed under the direction of Bear Braumoeller and with support from
## The Ohio State University's College of Social and Behavioral Sciences.
##
## boolean3 and is free software: you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation, either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program.  If not, see <http://www.gnu.org/licenses/>.
##
## ----------------------------------------------------------------------------


## Load default options.
.onLoad <- function(libname, pkgname) {
  boolean <- list(verbose = TRUE)
  options(boolean = boolean)
}


## Load message.
start_msg <- function(...) {
  lib   <- dirname(system.file(package = "boolean3"))
  vers  <- packageDescription("boolean3")$Version
  built <- packageDescription("boolean3")$Date
  
  packageStartupMessage(cat("## ", rep("-", 66), "\n", sep = ""))
  packageStartupMessage(cat("## boolean3: Modeling Causal Complexity"))
  packageStartupMessage(cat("## Version: ", vers))
  packageStartupMessage(cat("## Built:   ", built, "\n"))
  packageStartupMessage(cat("## When using this package, please cite:"))
  packageStartupMessage(cat("##    Braumoeller, Bear F. (2003) 'Causal Complexity and the Study"))
  packageStartupMessage(cat("##       of Politics'. Political Analysis 11(3): 209-233.\n"))
  packageStartupMessage(cat("## boolean3 was developed by Jason W. Morgan under the direction"))
  packageStartupMessage(cat("## of Bear Braumoeller with support from The Ohio State University's"))
  packageStartupMessage(cat("## College of Social and Behavioral Sciences. The package represents"))
  packageStartupMessage(cat("## a significant re-write of the original boolean implementation"))
  packageStartupMessage(cat("## developed by Bear Braumoeller, Ben Goodrich, and Jacob Kline."))
  packageStartupMessage(cat("## Please see the release notes and accompanying documentation for"))
  packageStartupMessage(cat("## details regarding changes made in this version.\n"))
  packageStartupMessage(cat("## ", rep("-", 66), sep = ""))
}

.onAttach <- function(...) {
  start_msg()
}
