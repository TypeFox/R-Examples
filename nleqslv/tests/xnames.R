
library(nleqslv)

f <- function(x) {
    y <-numeric(length(x))
    y[1] <- x[1]^2 + x[2]^3
    y[2] <- x[1] + 2*x[2] + 3
    y
}

# test named x-values
xstart <- c(a=1.0, b=0.5)
xstart

z <- nleqslv(xstart,f, control=list(trace=0))
all(names(z$x) == names(xstart))

# test named x-values
xstart <- c(u=1.0, 0.5)
xstart

z <- nleqslv(xstart,f, control=list(trace=0))
all(names(z$x) == names(xstart))
