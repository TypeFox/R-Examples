
library(nleqslv)

# Dennis Schnabel example 6.5.1 page 149
f <- function(x) {
    y <- numeric(2)
    y[1] <- x[1]^2 + x[2]^2 - 2
    y[2] <- exp(x[1]-1) + x[2]^3 - 2
    y
}

# check error handling in control argument
try(nleqslv(f,control=list(1e-3)))
try(nleqslv(f,control=list(f=1e-3)))
try(nleqslv(f,control=list(f=1e-7,b=1e-3)))
