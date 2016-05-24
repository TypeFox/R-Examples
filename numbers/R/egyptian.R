##
##  e g y p t i a n . R  Egyptian Fractions
##


egyptian_complete <- function(a, b, show = TRUE) {
    stopifnot(isNatural(a), isNatural(b), length(a) == 1, length(b) == 1)
	if (a >= b)
		stop("Rational number 'a/b' must be smaller than 1.")
	g <- GCD(a, b)
	if (g > 1) {
		warning("Arguments 'a', 'b' have a common divisor > 1.")
	    a <- a/g; b <- b/g
	}
	a0 <- a;
	b0 <- b

    noex <- 0  # no. of examples found

    # Define the first fraction: a/b - 1/x1 = a1/b1
	lb1 <- max(1 + div(b, a), 1)
	ub1 <- div(3*b, a)
	for (x1 in lb1:ub1) {
	    a1 <- a*x1 - b; b1 <- b*x1
		g <- GCD(a1, b1)
		a1 <- a1/g; b1 <- b1/g
		if (b1 == x1) next
		if (a1 == 1) {
            if (show)
			    cat("1/", x1, " + ", "1/", b1, "\n", sep = "")
            noex <- noex + 1
		} else {

            # Define the second fraction: a/ - 1/x1 - 1/x2 = a2/b2
            lb2 <- max(1 + div(b1, a1), x1)
			ub2 <- div(2*b1, a1)
			for (x2 in lb2:ub2) {
			    a2 <- a1*x2 - b1; b2 <- b1*x2
				g <- GCD(a2, b2)
				a2 <- a2/g; b2 <- b2/g
				if (x1 == x2 || b2 == x2) next
				if (a2 == 1) {
                    if (show)
				        cat("1/", x1, " + ", "1/", x2, " + ", "1/", b2, "\n", sep = "")
                    noex <- noex + 1
				}
			}
		}
	}
    invisible(noex)
}


egyptian_methods <- function(a, b) {
	stopifnot(isNatural(a), isNatural(b), length(a) == 1, length(b) == 1)
	if (a >= b)
		stop("Rational number 'a/b' must be smaller than 1.")
	g <- GCD(a, b)
	if (g > 1) {
		warning("Arguments 'a', 'b' have a common divisor > 1.")
	    a <- a/g; b <- b/g
	}

	a0 <- a; b0 <- b
	print_eg <- function(a, b, f) {
		cat(a, "/", b, " = ", sep="")
		cat("1/", f[1], sep="")
		for (i in 2:length(f)) {
			cat(" + 1/", f[i], sep = "")
		}
	}

	if (a == 1) {
		f <- c(b+1, b*(b+1))
		print_eg(a, b, f)
		cat("  (Egyptian fraction)\n")
		stop("Argument 'a' shall be an integer > 1.")
	}

	if (a == 2) {
	    if (isPrime(b) && b < 12) {
		    f <- c(b, b+1, b*(b+1))
		    print_eg(a, b, f)
    	    cat("  (Rhind Papyros)\n")
    	} else if (isPrime(b)) {
    	    f <- b * c(1, 2, 3, 6)
		    print_eg(a, b, f)
    	    cat("  (Rhind Papyros)\n")
		} else {
		    fctrs <- primeFactors(b)
		    if (length(fctrs) == 2 && all(isPrime(fctrs))) {
		        f <- (fctrs[1] + 1)/2 * c(1, fctrs[1]) * fctrs[2]
		        print_eg(a, b, f)
		        cat("  (Rhind Papyros)\n")
		    }
		}
	}

	# Fibonacci-Sylvester algorithm
	res <- c()
	while (a > 1) {
		s <- div(b, a)
		r <- mod(b, a)
		res <- c(res, s+1)
		a <- a - r
		b <- b*(s+1)
		g <- GCD(a, b)
		if (g != 1) {a <- a/g; b <- b/g}
	}
	res <- c(res, b)
	print_eg(a0, b0, res)
	cat("  (Fibonacci-Sylvester)\n")

	# Golomb-Farey algorithm
	a <- a0; b <- b0
	res <- c()
	while (a > 1) {
		p <- modinv(a, b)
		res <- c(res, p*b)
		a <- (p*a - 1)/b
		b <- p
		g <- GCD(a, b)
		if (g != 1) {a <- a/g; b <- b/g}
	}
	res <- rev(c(res, b))
	print_eg(a0, b0, res)
	cat("  (Golomb-Farey)\n")
}


##  TODO
# Extract the classical approaches into a subfunction called by main
# Add the following Papyros ideas
# (cf. http://en.wikipedia.org/wiki/Egyptian_fraction):
#
#   n/(p*q) = 1/(p*r) + 1/(q*r) with r = (p+q)/n, if possible (e.g., p = 1)
#   # example: 2/37 = 1/24 + 1/111 + 1/296, A = 24 (Ahmes)
#
#   2/p = 1/A + 1/(A*p)         where 2*A = p+1, p <= 20
#   # example: 2/13 = 1/7 + 1/91, A = 7
#
#   3/p = 1/p + 2/p, ...
#
# Add the following modern methods to egyptian_methods:
#
#   - greedy
#   - Bleicher-Erdos
#   - Tenenbaum-Yokota
#   - continued fractions
#
# define a 'best' measure like "geometricMean(f)", f the vector of 
# denominators, and introduce an option to retrieve this best solution.
#



