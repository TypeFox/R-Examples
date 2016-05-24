######################################################################
#
# zzz.R
#
# Written by Carter T. Butts <buttsc@uci.edu>.
#
# Last Modified 04/02/13
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/relevent package
#
# .onLoad is run when the package is loaded with library(relevent)
#
######################################################################

.onAttach <- function(libname, pkgname){
  temp<-packageDescription("relevent")
  msg<-paste(temp$Package,": ",temp$Title,"\n",
      "Version ",temp$Version,
      " created on ",
      temp$Date,".\n", sep="")
  msg<-paste(msg,"copyright (c) 2007, Carter T. Butts, University of California-
Irvine\n",sep="")
  msg<-paste(msg,'For citation information, type citation("relevent").\n')
  msg<-paste(msg,'Type help(package="relevent") to get started.\n')
  packageStartupMessage(msg)
}
