#  zzz.R
#
# This function is simply copied from mice package.

#------------------------------.onLoad-------------------------------
# .onLoad <- function(...){
#  d <- packageDescription("sirt")
#  cat("\n----------------------------\n")
#  packageStartupMessage(paste(d$Package," " , d$Version," (",d$Date,")",sep=""))
#  cat("See https://sites.google.com/site/alexanderrobitzsch/software\n")
#  cat("----------------------------\n")  
#  return()
# }
# on attach CDM
.onAttach <- function(libname,pkgname){
  d <- utils::packageDescription("sirt")
  d1 <- d$Version
# d1 <- "1.2"  
  nk <- paste( rep( " " , 20 - nchar(d1) ) , collapse="")
  packageStartupMessage("|---------------------------------------------------------",
		   "--------|\n"  ,
		paste("| " ,d$Package," " , d1 ," (",d$Date,")",sep="") , nk , 
		"                          |" , 
		"\n| Supplementary Item Response Theory                              |" ,
        "\n| Maintainer: Alexander Robitzsch <robitzsch@ipn.uni-kiel.de>     |" ,
		"\n| https://sites.google.com/site/alexanderrobitzsch/software       |",
		"\n|---------------------------------------------------" ,
		"--------------|\n" )
	}
	
	
version <- function(pkg="sirt"){
  lib <- dirname(system.file(package = pkg))
  d <- utils::packageDescription(pkg)
  return(paste(d$Package,d$Version,d$Date,lib))
}

# .First.lib <- function(lib, pkg){
#          library.dynam("sirt", package = pkg, lib.loc = lib)
#          return(invisible(0))
#        } 