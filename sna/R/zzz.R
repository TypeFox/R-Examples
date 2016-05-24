######################################################################
#
# zzz.R
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 2/28/13
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/sna package
#
# .onLoad is run when the package is loaded with library(sna)
#
######################################################################

.onLoad <- function(libname, pkgname){
  library.dynam("sna", package=pkgname, lib.loc=libname)
}

.onAttach <- function(libname, pkgname){
  temp<-packageDescription("sna")
  msg<-paste(temp$Package,": ",temp$Title,"\n",
      "Version ",temp$Version,
      " created on ",
      temp$Date,".\n", sep="")
  msg<-paste(msg,"copyright (c) 2005, Carter T. Butts, University of California-Irvine\n",sep="")
  msg<-paste(msg,'For citation information, type citation("sna").\n')
  msg<-paste(msg,'Type help(package="sna") to get started.\n')
  packageStartupMessage(msg)
}

