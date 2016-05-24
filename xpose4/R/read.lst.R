# Xpose 4
# An R-based population pharmacokinetic/
# pharmacodynamic model building aid for NONMEM.
# Copyright (C) 1998-2004 E. Niclas Jonsson and Mats Karlsson.
# Copyright (C) 2005-2008 Andrew C. Hooker, Justin J. Wilkins, 
# Mats O. Karlsson and E. Niclas Jonsson.
# Copyright (C) 2009-2010 Andrew C. Hooker, Mats O. Karlsson and 
# E. Niclas Jonsson.

# This file is a part of Xpose 4.
# Xpose 4 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation, either version 3
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program.  A copy can be cound in the R installation
# directory under \share\licenses. If not, see http://www.gnu.org/licenses/.

read.lst <- function(filename) {

  ## The function determines the NONMEM version used to produce the
  ## list file and invokes the appropriate read.lst function.

  listfile <- scan(filename, sep = "\n", what = character(),quiet=TRUE)

  ## Match the VERSION string
  versionline <- grep("1NONLINEAR", listfile)

  if(is.null(version$language)){
    cat("need to use R for this version of Xpose")
    ##&&
    ##  platform() == "WIN386" &&
    ##  version$major <6) {
    ## versionVIpat  <- "*VERSION*VI\\s*"        #Not tested
    ## versionVIIpat <- "*VERSION*7.*"       #Not tested
  } else {
    versionVIpat  <- "VERSION VI\\s"
    versionVIIpat <- "VERSION 7."
  }

  versionVIpatline  <- grep(versionVIpat, listfile)
  versionVIIpatline <- grep(versionVIIpat, listfile)

  ## Check that we found a NONMEM version
  if(length(versionVIpatline) == 0 && length(versionVIIpatline)==0) stop("Can not establish NONMEM version\n")

  NMversion <- NA
  ## Check which version we are dealing with
  if(length(versionVIpatline) !=0 ) {
    if(any(versionline == versionVIpatline)) NMversion <- 6
  }

  if(length(versionVIIpatline) !=0 && is.na(NMversion)) {
    if(any(versionline == versionVIIpatline)) NMversion <- 7
  }

  if(is.na(NMversion)) stop("Can not establish NONMEM version\n")

  if(NMversion==6) {
    return(read.lst6(filename))
  } else {
    return(read.lst7(filename))
  }


}
