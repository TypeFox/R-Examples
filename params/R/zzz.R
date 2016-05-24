.onAttach <- function(lib, pkg){
}

.onDeatch <- function(lib, pkg) {
	#ops <- options()
	#ops <- ops[grep("SaturnV_", names(ops))]
	#options(ops)
}


.onLoad <- function(lib, pkg){
	fls = system.file("conf/params.conf", package = "params")
	suppressMessages(load_opts(fls, check = FALSE))
}
