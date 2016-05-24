.onLoad <-function(libname, pkgname){
  library.dynam("networksis", package=pkgname, lib.loc=libname)
}

.onAttach <- function(libname, pkgname){
  temp<-packageDescription("networksis")
  msg<-paste(temp$Package,": ",temp$Title,"\n",
      "Version ",temp$Version,
      " created on ",
      temp$Date,".\n", sep="")
  msg<-paste(msg,"copyright (c) 2008, Ryan Admiraal, Murdoch University\n",sep="")
  msg<-paste(msg,'Based on "statnet" project software (statnet.org).\n')
  msg<-paste(msg,'For license and citation information see statnet.org/attribution\n')
  msg<-paste(msg,'For citation information, type citation("networksis").\n')
  msg<-paste(msg,'Type help("networksis-package") to get started.\n')
  packageStartupMessage(msg)
}
