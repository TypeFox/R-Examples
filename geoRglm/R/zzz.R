".onAttach" <- function(lib, pkg)
{
  messages.screen <- ifelse(is.null(getOption("geoR.messages")), TRUE, getOption("geoR.messages"))
  packageStartupMessage("---------------------------------------------------------\n")
  pkg.info <- drop(read.dcf(file=system.file("DESCRIPTION", package="geoRglm"), fields=c("Title","Version","Date")))
  packageStartupMessage(pkg.info["Title"])
  packageStartupMessage("\n")
  packageStartupMessage(paste("geoRglm version ", pkg.info["Version"], " (", pkg.info["Date"], ") is now loaded\n", sep=""))
  packageStartupMessage("-----------------------------------------------------------\n")
  packageStartupMessage("\n")
  return(invisible(0))
}


## geoR functions, that are now copied into geoRglm
##".bilinearformXAY"<- function(...) geoR:::.bilinearformXAY(...)
##".check.locations"<- function(...) geoR:::.check.locations(...)
##".cond.sim"<- function(...) geoR:::.cond.sim(...)
##".cor.number"<- function(...) geoR:::.cor.number(...)
##".diagquadraticformXAX"<- function(...) geoR:::.diagquadraticformXAX(...)
##".geoR_inout"<- function(...) geoR:::.geoR_inout(...)
##".ldots.set"<- function(...) geoR:::.ldots.set(...)
##".solve.geoR"<- function(...) geoR:::.solve.geoR(...)
