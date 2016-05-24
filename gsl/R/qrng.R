qrng_alloc <- function(type = c("niederreiter_2", "sobol"), dim) {
    type <- switch(match.arg(type),
    		niederreiter_2 = 0,
    		sobol = 1)
    .Call("qrng_alloc", type, dim, PACKAGE = "gsl")
}

qrng_clone <- function(q) .Call("qrng_clone", q, PACKAGE = "gsl")

qrng_init <- function(q) .Call("qrng_init", q, PACKAGE = "gsl")

qrng_name <- function(q) .Call("qrng_name", q, PACKAGE = "gsl")

qrng_size <- function(q) .Call("qrng_size", q, PACKAGE = "gsl")

qrng_get <- function(q, n = 1) 
    matrix(.Call("get_n", q, n, PACKAGE = "gsl"), nrow = n, byrow = TRUE)
    
