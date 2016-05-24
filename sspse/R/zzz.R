######################################################################
# copyright (c) 2009, Krista J. Gile, University of Massachusetts - Amherst
#                     Mark S. Handcock, University of California - Los Angeles
# 
# For license and citation information see
#    http://statnet.org/attribution
#
# We have invested a lot of time and effort in creating 'statnet',
# for use by other researchers. We require that the attributions
# in the software are retained (even if only pieces of it are used),
# and that there is attribution when the package is loaded (e.g., via
# "library" or "require"). 
######################################################################
# File name: zzz.R
######################################################################
#
# .First.lib is run when the package is loaded.
#
######################################################################

.onLoad <-function(libname, pkgname){
  library.dynam("sspse", package=pkgname, lib.loc=libname)
}

.onAttach <- function(libname, pkgname){
  temp<-packageDescription("sspse")
  msg<-paste(temp$Package,": ",temp$Title,"\n",
      "Version ",temp$Version,
      " created on ",
      temp$Date,".\n", sep="")
  msg<-paste(msg,"copyright (c) 2014, Krista J. Gile, University of Massachusetts - Amherst\n",
"                    Mark S. Handcock, University of California - Los Angeles\n",sep="")
  msg<-paste(msg,'For citation information, type citation("sspse").\n')
  msg<-paste(msg,'Type help("sspse-package") to get started.\n')
  packageStartupMessage(msg)
}
