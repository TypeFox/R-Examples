# R. Baker Kearfott, Some tests of Generalized Bisection,
# ACM Transactions on Methematical Software, Vol. 13, No. 3, 1987, pp 197-220

# A high-degree polynomial system (section 4.3 Problem 12)
# There are 12 real roots (and 126 complex roots to this system!)

library(nleqslv)

hdp <- function(x) {
    f <- numeric(length(x))
    f[1] <- 5 * x[1]^9 - 6 * x[1]^5 * x[2]^2 + x[1] * x[2]^4 + 2 * x[1] * x[3]
    f[2] <- -2 * x[1]^6 * x[2] + 2 * x[1]^2 * x[2]^3 + 2 * x[2] * x[3]
    f[3] <- x[1]^2 + x[2]^2 - 0.265625
    f
}


N <- 40
set.seed(123)
xstart <- matrix(runif(3*N,min=-1,max=1), N, 3)  # N starting values, each of length 3

ans <- searchZeros(xstart,hdp, method="Broyden",global="dbldog")
nrow(ans$x) == 12
all(ans$xfnorm <= 1e-10)

zans <- searchZeros(ans$xstart,hdp, method="Broyden",global="dbldog")
length(zans$idxcvg) == 12
