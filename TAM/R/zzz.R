#  zzz.R
#
# This function is simply copied from mice package.

#------------------------------.onLoad-------------------------------
#.onLoad <- function(...){
#  d <- packageDescription("TAM")
#  cat("\n............................\n")
#  packageStartupMessage(paste(d$Package," " , d$Version," (",d$Date,")",sep=""))
#  cat("............................\n")  
#  return()
#}
version <- function(pkg="TAM"){
  lib <- dirname(system.file(package = pkg))
  d <- utils::packageDescription(pkg)
  return(paste(d$Package,d$Version,d$Date,lib))
}
# on attach TAM
.onAttach <- function(libname,pkgname){
  d <- utils::packageDescription("TAM")
  base::packageStartupMessage(  # "::...........................::\n",
		paste0("* " , d$Package," " , d$Version," (",d$Date,")" )  # , 
#          paste0(rep(" ", max(8-nchar(d$Version),0)), collapse=""), " ::",sep="") ,
#		paste("\n:: Test Analysis Modules     ::") ,
#		"\n::...........................::\n" 
						)
}

#*************************************************
# R session info
Rsessinfo <- function(){
    si <- Sys.info()
    si2 <- utils::sessionInfo()
    paste0( si2$R.version$version.string , " " , si2$R.version$system 
             , " | nodename = " , si["nodename"] , " | login = " , si["login"] )
            }