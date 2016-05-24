######################################################################
#
# zzz.R
#
# Written by Christopher Steven Marcum <cmarcum@uci.edu>.
# Adapted from Carter T. Butts
#
# Last Modified 04/02/13
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/relevent package
#
# .onLoad is run when the package is loaded with library(informR)
#
######################################################################

.onAttach <- function(libname, pkgname){
  temp<-packageDescription("informR")
  msg<-paste(temp$Package,": ",temp$Title,"\n",
      "Version ",temp$Version,
      " created on ",
      temp$Date,".\n", sep="")
  msg<-paste(msg,"copyright (c) 2010, Christopher Steven Marcum, University of California-
Irvine\n",sep="")
  msg<-paste(msg,'For citation information, type citation("informR").\n')
  msg<-paste(msg,'Type help(package="informR") to get started.\n')
  packageStartupMessage(msg)
}
