.onAttach <- function(libname, pkgname){
  temp<-packageDescription("reldist")
  msg<-paste(temp$Package,": ",temp$Title,"\n",
      "Version ",temp$Version,
      " created on ",
      temp$Date,".\n", sep="")
  msg<-paste(msg,"copyright (c) 2003, Mark S. Handcock, University of California-Los Angeles\n",sep="")
  msg<-paste(msg,'For citation information, type citation("reldist").\n')
# msg<-paste(msg,'Type help("reldist-package") to get started.\n')
  msg<-paste(msg,'Type help(package="reldist") to get started.\n')
  packageStartupMessage(msg)
}
