## welcome message

.onAttach <- function(libname, pkgname){
	suppressWarnings(descr <- utils::packageDescription("PST"))

	builtDate <- strsplit(strsplit(descr$Built, ";")[[1]][3], " ")[[1]][2]

	packageStartupMessage("\n",descr$Package," version ", descr$Version, " (Built: ", builtDate, ")")
	packageStartupMessage("Website: ", descr$URL)
	## packageStartupMessage("Please type 'citation(\"PST\")' for citation information.\n")
}
