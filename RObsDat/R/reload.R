reload <- function(){
	try(detach("package:RObsDat", unload=TRUE))
	ll = getOption("lib.loc", default=NULL)
	ll.string <- ""
	if(!is.null(ll)) ll.string <- paste("-l", ll)
	system(paste("R CMD INSTALL",ll.string," RObsDat"))
	library(RObsDat, lib.loc=ll)
}
