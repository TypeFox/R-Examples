"rng_alloc" <- function(type = 
			c("mt19937",
			"ranlxs0",
			"ranlxs1",
			"ranlxs2",
			"ranlxd1",
			"ranlxd2",
			"ranlux",
			"ranlux389",
			"cmrg",
			"mrg",
			"taus",
			"taus2",
			"gfsr4",
			"minstd")) {

	type <- switch( match.arg(type),
			mt19937 = 0,
			ranlxs0 = 1,
			ranlxs1 = 2,
			ranlxs2 = 3,
			ranlxd1 = 4,
			ranlxd2 = 5,
			ranlux = 6,
			ranlux389 = 7,
			cmrg = 8,
			mrg = 9,
			taus = 10,
			taus2 = 11,
			gfsr4 = 12,
			minstd = 13
					)
	retval <- .Call("rng_alloc", type, PACKAGE = "gsl")
	class(retval) <- c("gsl_rng", class(retval))
	retval
}

"rng_set" <-
function(r, seed)
{
    .Call("rng_set", r, seed, PACKAGE = "gsl")
}
"rng_clone" <- function(r) {
	retval <- .Call("rng_clone", r, PACKAGE = "gsl")
	class(retval) <- "gsl_rng"
	retval
}
"rng_name" <-
function(r)
{
    .Call("rng_name", r, PACKAGE = "gsl")
}
"rng_min" <-
function(r)
{
    .Call("rng_min", r, PACKAGE = "gsl")
}
"rng_max" <-
function(r)
{
    .Call("rng_max", r, PACKAGE = "gsl")
}
"rng_get" <-
function(r, length)
{
    .Call("rng_get", r, length, PACKAGE = "gsl")
}
"rng_uniform" <-
function(r, length)
{
    .Call("rng_uniform", r, length, PACKAGE = "gsl")
}
"rng_uniform_pos" <-
function(r, length)
{
    .Call("rng_uniform_pos", r, length, PACKAGE = "gsl")
}
"rng_uniform_int" <- function(r, N, length)
{
    if( ! (N > 0) )
      stop("N needs to be positive")
    .Call("rng_uniform_int", r, N, length, PACKAGE = "gsl")
}
