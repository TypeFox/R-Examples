###
### $Id: factors.R 51 2014-02-05 21:22:28Z plroebuck $
###


##-----------------------------------------------------------------------------
test.factors <- function(input, expected) {
    output <- do.call(getFromNamespace("factors", "matlab"), input)
    identical(output, expected)
}

factors.expected.n2  <- 2
factors.expected.n3  <- 3
factors.expected.prm <- 999983
factors.expected.pr2 <- c(9999889, 9999901)
factors.expected.prp <- c(65003, 65003)
factors.expected.nn  <- c(2, 2, 2, 2, 2, 3, 3, 3, 3, 5, 5, 5, 7, 7, 11)
factors.expected.nm  <- c(99989, 99991, 100003)
factors.expected.n32 <- c(3, 5, 17, 257, 65537)

test.factors(list(n=2), factors.expected.n2)
test.factors(list(n=3), factors.expected.n3)
test.factors(list(n=999983), factors.expected.prm)
#test.factors(list(n=9999889*9999901), factors.expected.pr2)
test.factors(list(n=4225390009), factors.expected.prp)
test.factors(list(n=2^5 * 3^4 * 5^3 * 7^2 * 11), factors.expected.nn)
test.factors(list(n=99989*99991*100003), factors.expected.nm)
test.factors(list(n=2^32-1), factors.expected.n32)

