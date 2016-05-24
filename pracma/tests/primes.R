###
### PRIMES.R  +++ Test suite +++
###


test.primes <- function(input, expected) {
    output <- do.call(getFromNamespace("primes", "pracma"), input)
    identical(output, expected)
}

primes.expected.n1  <- NULL
primes.expected.n2  <- 2
primes.expected.n3  <- c(2, 3)
primes.expected.n5  <- c(2, 3, 5)
primes.expected.n8  <- c(2, 3, 5, 7)
primes.expected.101 <- c( 2,  3,  5,  7, 11, 13, 17, 19, 23, 29, 31, 37, 41,
                         43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,101)

test.primes(list(n=1), primes.expected.n1)
test.primes(list(n=2), primes.expected.n2)
test.primes(list(n=3), primes.expected.n3)
test.primes(list(n=5), primes.expected.n5)
test.primes(list(n=8), primes.expected.n8)
test.primes(list(n=101), primes.expected.101)
