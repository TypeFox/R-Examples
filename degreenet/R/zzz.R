.onLoad <-function(libname, pkgname){
  library.dynam("degreenet", package=pkgname, lib.loc=libname)
}

.onAttach <- function(libname, pkgname){
  temp<-packageDescription("degreenet")
  msg<-paste(temp$Package,": ",temp$Title,"\n",
      "Version ",temp$Version,
      " created on ",
      temp$Date,".\n", sep="")
  msg<-paste(msg,"copyright (c) 2013, Mark S. Handcock, University of California - Los Angeles\n",sep="")
  msg<-paste(msg,'Based on "statnet" project software (statnet.org).\n')
  msg<-paste(msg,'For license and citation information see statnet.org/attribution\n')
  msg<-paste(msg,'For citation information, type citation("degreenet").\n')
  msg<-paste(msg,'Type help("degreenet-package") to get started.\n')
  packageStartupMessage(msg)
}
