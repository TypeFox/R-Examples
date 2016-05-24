.onLoad <- function(lib, pkg)
{
	library.dynam("rngwell19937", pkg, lib)
	RNGkind("user-supplied")
}

.onUnload <- function(libpath) {
	RNGkind("default")
	library.dynam.unload("rngwell19937", libpath)
}

set.resolution <- function(resolution)
{
	if (resolution == 53 || resolution == 32) {
		aux <- .C("set_resolution",
			as.logical(resolution == 53),
			PACKAGE="rngwell19937")
	} else {
		cat("supported resolutions are 53 and 32 bits\n")
	}
	invisible(NULL)
}

set.initialization <- function(initialization)
{
	if (initialization == "mrg32k5a" || initialization == "sfmt") {
		aux <- .C("set_initialization",
			as.integer(initialization == "mrg32k5a"),
			PACKAGE="rngwell19937")
	} else {
		cat("supported initializations are \"mrg32k5a\" and \"sfmt\"\n")
	}
	invisible(NULL)
}

set.vector.seed <- function(seed)
{
	stopifnot(all(seed >= 0))
	stopifnot(all(seed < 2^32))
	stopifnot(all(floor(seed) == seed))
	RNGkind("user-supplied")
	tmp <- .C("init_vector_mrg32k5a",
		length(seed),
		as.double(seed),
		integer(length(seed)),
		new.state = integer(625),
		PACKAGE="rngwell19937")
	.Random.seed[2:626] <<- tmp$new.state
	invisible(NULL)
}

