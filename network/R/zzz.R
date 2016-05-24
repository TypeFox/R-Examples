######################################################################
#
# zzz.R
#
# Written by Carter T. Butts <buttsc@uci.edu>.
#
# Last Modified 03/01/12
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/network package
#
# .onLoad is run when the package is loaded with library(network)
#
######################################################################

.onLoad <-function(libname, pkgname){
  library.dynam("network", package=pkgname, lib.loc=libname)
}

.onAttach <- function(libname, pkgname){
  temp<-packageDescription("network")
  msg<-paste(temp$Package,": ",temp$Title,"\n",
      "Version ",temp$Version,
      " created on ",
      temp$Date,".\n", sep="")
  msg<-paste(msg,"copyright (c) 2005, Carter T. Butts, University of California-Irvine\n",
"                    Mark S. Handcock, University of California -- Los Angeles\n",
"                    David R. Hunter, Penn State University\n",
"                    Martina Morris, University of Washington\n",
"                    Skye Bender-deMoll, University of Washington\n",sep="")
  msg<-paste(msg,'For citation information, type citation("network").\n')
  msg<-paste(msg,'Type help("network-package") to get started.\n')
  packageStartupMessage(msg)
}