library(randtoolbox)

x <- halton(1, 205)
x <- halton(1, 5)
print(halton(1, 5))

options(digits=15)
cbind( 1/get.primes(500), as.vector(halton(1, 500) ) )

