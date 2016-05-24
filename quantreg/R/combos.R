"combos" <- function(n,p){
    if(length(n) != 1){ 
        n <- n[1]
        warning("Using first element as n")
        }
    if(length(p) != 1){ 
        p <- p[1]
        warning("Using first element as p")
        }
    if(n != as.integer(n)){
        warning("Coercing n to integer")
        n <- as.integer(n)
        }
    if(p != as.integer(p)){
        warning("Coercing p to integer")
        p <- as.integer(p)
        }
    if(p > n) stop("p is greater than n")
    m <- choose(n,p)
    z <- .Fortran("combin",
		as.integer(n),
		as.integer(p),
		as.integer(m),
		a = integer(p*m),
		integer(n),
		integer(n),
		integer(n),
		PACKAGE = "quantreg")
	matrix(z$a,p)
	}
