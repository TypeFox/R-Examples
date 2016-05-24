.onAttach <- function(lib, pkg) {
#	library.dynam("spgwr", pkg, lib)
	packageStartupMessage(paste(
            "NOTE: This package does not constitute approval of GWR",
            "as a method of spatial analysis; see example(gwr)\n", sep="\n"),
            appendLF = FALSE)
}
