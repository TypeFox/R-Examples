
.onAttach <- function(libname, pkgname){
	suppressWarnings(descr <- packageDescription("TraMineRextras"))
	if(packageVersion("TraMineRextras")$minor %% 2 == 0) {
		state <- "stable"
	}
	else {
		state <- "development"
	}
	builtDate <- strsplit(strsplit(descr$Built, ";")[[1]][3], " ")[[1]][2]
	packageStartupMessage("TraMineRextras ", state, " version ", descr$Version, " (Built: ", builtDate, ")")
	packageStartupMessage("Functions provided by this package are still in test")
	packageStartupMessage("    and subject to changes in future releases.")
}
