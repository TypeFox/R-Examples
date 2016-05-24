
.onLoad <- function(lib, pkg) {
	library.dynam("rstream", pkg, lib)
	.rstream.init()
	.rstream.lecuyer.init()
	.rstream.mrg32k3a.init()
	.rstream.runif.init()
}

.onUnload <- function(libpath) {
	library.dynam.unload("rstream", libpath)
}
