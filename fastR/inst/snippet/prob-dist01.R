require(MASS)                        # for fractions()
kernel <- function(x) { (x-2)*(x+2) * as.numeric(x >=-2 & x <= 2) }
k <- 1 / integrate(kernel,-2,2)$value; k
f <- function(x) { k * kernel(x) }
fractions(k)
integrate(f,-2,2)                    # check that we have pdf
integrate(f,0,2)
fractions(integrate(f,0,2)$value)
integrate(f,1,2)
fractions(integrate(f,1,2)$value)
integrate(f,-1,1)
fractions(integrate(f,-1,1)$value)
