## .First.lib <- function(lib, pkg) {library.dynam("TraMineR", pkg, lib)}

.onAttach <- function(libname, pkgname){
	suppressWarnings(descr <- utils::packageDescription("TraMineR"))
	Tver <- extract.ver(descr$Version)
	if(as.numeric(Tver[2])%%2==0) {
		state <- "stable"
	}
	else {
		state <- "development"
	}
	builtDate <- strsplit(strsplit(descr$Built, ";")[[1]][3], " ")[[1]][2]
	packageStartupMessage("\n",descr$Package," ", state, " version ", descr$Version, " (Built: ", builtDate, ")")
	packageStartupMessage("Website: ", descr$URL)
	packageStartupMessage("Please type 'citation(\"TraMineR\")' for citation information.\n")
}
